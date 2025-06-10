#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define LATGP latgp
#define LATG2P latg2p
#define LONF lonf
      subroutine mpgp2f(a,lonf2p,latg2p,b,lonf2,latg2,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mpgp2f
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  transpose (ip,jp,kf) to (if,jf,kf)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpgp2f(a,lonf2p,latg2p,b,lonf2,latg2,ntotal)
c
c    input argument lists:
c	b	- real (lonf2p,latg2p,ntotal) partial field for each pe 
c	lonf2p	- integer partial grid in longitude time 2
c	latg2p	- integer partial grid in latitude divide 2
c	ntotal	- integer total set of fields
c
c    output argument list:
c	a	- real (lonf2,latg2,ntotal) full field 
c	lonf2	- integer total grid in longitude time 2
c	latg2	- integer total grid in latitude divide 2
c 
c subprograms called:
c   mpi_gatherv	- gather message from all pe to master
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
      integer lonf2p,latg2p,lonf2,latg2,ntotal
     1       ,ii,jj,n,i,j,k,mk,ierr,lonlend,lonstrd,lonff
      real    a,b
      dimension a(lonf2p,latg2p,ntotal),b(lonf2,latg2,ntotal)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(lonf2p*latg2p*ntotal))
      allocate(tmprcv(lonf2 *latg2 *ntotal))
      allocate(len(0:npes_-1))
      allocate(loc(0:npes_-1))
c
c
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
#ifdef REDUCE_GRID
          jj=0
          do j=1,latlen(n)
            jj=jj+lonfdp(j,n)
          enddo
          len(n)=ntotal*jj*2
#else
          len(n)=ntotal*lonlen(n)*2*latlen(n)
#endif
          mk=loc(n)+len(n)
        enddo
      endif
c
#ifndef REDUCE_GRID
      lonlend=lonlen(mype)
#endif
      mk=0
      do k=1,ntotal
        do j=1,latlen(mype)
#ifdef REDUCE_GRID
          lonlend=lonfdp(j,mype)
#endif
          do i=1,lonlend
            mk=mk+1
            tmpsnd(mk)=a(i,j,k)
            mk=mk+1
            tmpsnd(mk)=a(i+lonlend,j,k)
          enddo
        enddo
      enddo
      len(mype)=mk
c
      call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,
     1     tmprcv,len(0),loc(0),MPI_REAL,
     2     0,MPI_COMM_GSM,ierr)
c
      if( mype.eq.0 ) then
        mk=0
        do n=0,npes-1
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
                b(ii,jj,k)=tmprcv(mk)
                ii=ii+lonff
                mk=mk+1
                b(ii,jj,k)=tmprcv(mk)
              enddo
            enddo
          enddo
        enddo
      endif
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
