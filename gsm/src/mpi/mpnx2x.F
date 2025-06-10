#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATGP latgp
#define LATG2P latg2p
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
      subroutine mpnx2x(a,lonf2p,lota,b,lonf2,lotb,latg2p,
     &                  leva,levb,nvar)
c$$$  subprogram documentation block
c
c subprogram: 	mpnx2x
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:   transpose (ip,k,jp) to (if,k,jp)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpnx2x(a,lonf2p,lota,b,lonf2,lotb,latg2p,
c        &                  leva,levb,nvar)
c
c    input argument lists:
c	a	- real (lonf2p,lota,latg2p) partial field in i j
c	lonf2p	- integer partial longitude grid x 2
c	latg2p	- integer partial latitude grid / 2
c	lonf2	- integer total longitude grid x 2
c	leva	- integer starting vertical layer for a
c	levb	- integer starting vertical layer for b
c	nvar	- integer total set of fields
c
c    output argument list:
c	b	- real (lonf2,lotb,latg2p) partial field in j k
c 
c subprograms called:
c   mpi_alltoallv  - send and receive from all pe in the same comm
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
      integer lonf2p,lonf2,latg2p,lota,lotb,nvar,leva,levb,status
     1       ,i,j,k,ii,ipe,ipe0,ioff,noffa,noffb,n,mn,len,ierr
     1       ,lonf2d,lonfp,lonff,latlon
      real    a,b
      dimension a(lonf2p,lota,latg2p),b(lonf2,lotb,latg2p)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::lensnd(:),lenrcv(:)
      integer,allocatable::locsnd(:),locrcv(:)
c
c option for 1-d decomposition
      if( ncol.eq.1 ) then
        do j=1,latg2p
          do k=1,lota
            do i=1,lonf2
              b(i,k,j) = a(i,k,j)
            enddo
          enddo
        enddo
        return
      endif
c
c otherwise
      allocate(tmpsnd(lonf2p*latg2p*nvar*ncol))
      allocate(tmprcv(lonf2p*latg2p*nvar*ncol))
      allocate(lensnd(ncol))
      allocate(lenrcv(ncol))
      allocate(locsnd(ncol))
      allocate(locrcv(ncol))
c
      ipe0=int(mype/ncol)*ncol
c
      mn=0
#ifndef REDUCE_GRID
      lonfp=lonlen(mype)
#endif
      do ii=1,ncol
        locsnd(ii)=mn
        do n=1,nvar
          noffa=leva-1+n
          do j=1,latlen(mype)
#ifdef REDUCE_GRID
            lonfp=lonfdp(j,mype)
#endif
            do i=1,lonfp
              mn=mn+1
              tmpsnd(mn)=a(i      ,noffa,j)
              mn=mn+1
              tmpsnd(mn)=a(i+lonfp,noffa,j)
            enddo
          enddo
        enddo
        lensnd(ii)=mn-locsnd(ii)
      enddo
c
      mn=0
      do ii=1,ncol
        locrcv(ii)=mn
        ipe=(ii-1)+ipe0
#ifdef REDUCE_GRID
        latlon=0
        do j=1,latlen(mype)
          latlon=latlon+lonfdp(j,ipe)
        enddo
#else
        latlon=latlen(mype)*lonlen(ipe)
#endif
        lenrcv(ii)=latlon*2*nvar
        mn=mn+lenrcv(ii)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   comm_row,ierr)
c
      mn=0
      do ii=1,ncol
        ipe=(ii-1)+ipe0
#ifndef REDUCE_GRID
        lonff=lonf_
        lonfp=lonlen(ipe)
        ioff=lonstr(ipe)-1
#endif
        do n=1,nvar
          noffb=levb-1+n
          do j=1,latlen(mype)
#ifdef REDUCE_GRID
            lonff=lonfd(latdef(latstr(mype)+j-1))
            lonfp=lonfdp(j,ipe)
            ioff=lonfds(j,ipe)-1
#endif
            do i=1,lonfp
              mn=mn+1
              b(i+ioff      ,noffb,j)=tmprcv(mn)
              mn=mn+1
              b(i+ioff+lonff,noffb,j)=tmprcv(mn)
            enddo
          enddo
        enddo
      enddo
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(lensnd)
      deallocate(lenrcv)
      deallocate(locsnd)
      deallocate(locrcv)
c
      return
      end
