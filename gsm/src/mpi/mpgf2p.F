#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATGP latgp
#define LATG2P latg2p
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
      subroutine mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mpgf2p
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (if,jf,kf) to (ip,jp,kf)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
c
c    input argument lists:
c	a	- real (lonf2,latg2,ntotal) full field 
c	lonf2	- integer total grid in longitude time 2
c	latg2	- integer total grid in latitude divide 2
c	ntotal	- integer total set of fields
c
c    output argument list:
c	b	- real (lonf2p,latg2p,ntotal) partial field for each pe 
c	lonf2p	- integer partial grid in longitude time 2
c	latg2p	- integer partial grid in latitude divide 2
c 
c subprograms called:
c   mpi_scatterv	- scatter message from master pe to all pe
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
c
      implicit none
#include <paramodel.h>
#include <commpi.h>
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
      integer lonf2,latg2,lonf2p,latg2p,ntotal
     1       ,ii,jj,n,mk,i,j,k,ierr,lonlend,lonstrd,lonff
      real    a,b
      dimension a(lonf2,latg2,ntotal),b(lonf2p,latg2p,ntotal)
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(lonf2 *latg2 *ntotal))
      allocate(tmprcv(lonf2p*latg2p*ntotal))
      allocate(len(0:npes_-1))
      allocate(loc(0:npes_-1))
c
c
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
#ifndef REDUCE_GRID
          lonlend=lonlen(n)
          lonstrd=lonstr(n)-1
          lonff=lonf_
#endif
          do k=1,ntotal
            do j=1,latlen(n)
              jj=j+latstr(n)-1
              jj=latdef(jj)
#ifdef REDUCE_GRID
              lonlend=lonfdp(j,n)
              lonstrd=lonfds(j,n)-1
              lonff=lonfd(jj)
#endif
              do i=1,lonlend
                ii=i+lonstrd
                mk=mk+1
                tmpsnd(mk)=a(ii,jj,k)
                ii=ii+lonff
                mk=mk+1
                tmpsnd(mk)=a(ii,jj,k)
              enddo
            enddo
          enddo
          len(n)=mk-loc(n)
        enddo
      else
#ifdef REDUCE_GRID
        mk=0
        do j=1,latlen(mype)
          mk=mk+lonfdp(j,mype)
        enddo
        len(mype)=ntotal*mk*2
#else
        len(mype)=ntotal*lonlen(mype)*2*latlen(mype)
#endif
      endif
c
      call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,
     1     tmprcv,len(mype),MPI_REAL,0,MPI_COMM_GSM,ierr)
c
      mk=0
#ifndef REDUCE_GRID
      lonlend=lonlen(mype)
#endif
      do k=1,ntotal
        do j=1,latlen(mype)
#ifdef REDUCE_GRID
          lonlend=lonfdp(j,mype)
#endif
          do i=1,lonlend
            mk=mk+1
            b(i        ,j,k)=tmprcv(mk)
            mk=mk+1
            b(i+lonlend,j,k)=tmprcv(mk)
          enddo
          do i=lonlend*2+1,lonf2p
            b(i,j,k)=tmprcv(mk)
          enddo
        enddo
        do j=latlen(mype)+1,latg2p
          do i=1,lonf2p
            b(i,j,k)=tmprcv(mk)
          enddo
        enddo
      enddo
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
