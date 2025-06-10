#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATGP latgp
#define LATG2P latg2p
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
      subroutine mpnk2nx(a,lonf2,lota,
     &                   b,lonf2p,lotb,latg2p,levsp,levs,
     &                   leva,levb,nvar)
c$$$  subprogram documentation block
c
c subprogram: 	mpnk2nx
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  transpose (if,kp) to (ip,kf)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpnk2nx(a,lonf2,lota,
c        &                   b,lonf2p,lotb,latg2p,levsp,levs,
c        &                   leva,levb,nvar)
c
c    input argument lists:
c	a	- real (lonf2,lota,latg2p) partial field in j k
c	lonf2	- integer total longitude grid x 2
c	lota	- integer total field of nvar*levsp
c	lonf2p	- integer partial longitude grid x 2
c	lotb	- integer total field of nvar*levs
c	latg2p	- integer partial latitude grid / 2
c	levsp	- integer partial vertical grid 
c	levs	- integer total vertical grid 
c	leva	- integer start vertical grid for a
c	levb	- integer start vertical grid for b
c	nvar	- integer total set of fields
c
c    output argument list:
c	b	- real (lonf2p,lotb,latg2p) partial field in i j
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
      integer lonf2,levsp,latg2p,lonf2p,levs,nvar,status
     1       ,lota,lotb,leva,levb,ioff,koff,lonff,lonf2d
     1       ,i,j,k,ii,ipe,ipe0,n,m,mn,len,kp,kf,ierr
     1       ,latlon,lonfp
      real    a,b
      dimension a(lonf2,lota,latg2p),b(lonf2p,lotb,latg2p)
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
      allocate(tmpsnd(lonf2p*latg2p*levsp*nvar*ncol))
      allocate(tmprcv(lonf2p*latg2p*levsp*nvar*ncol))
      allocate(lensnd(ncol))
      allocate(lenrcv(ncol))
      allocate(locsnd(ncol))
      allocate(locrcv(ncol))
c
c
      ipe0=int(mype/ncol)*ncol
c cut in x
      mn=0
      do ii=1,ncol
        locsnd(ii)=mn
        ipe=(ii-1)+ipe0
#ifndef REDUCE_GRID
        lonff=lonf_
        lonfp=lonlen(ipe)
        ioff=lonstr(ipe)-1
#endif
        do n=1,nvar
          kp=(n-1)*levsp+leva-1
          do j=1,latlen(mype)
#ifdef REDUCE_GRID
            lonff=lonfd(latdef(latstr(mype)+j-1))
            lonfp=lonfdp(j,ipe)
            ioff=lonfds(j,ipe)-1
#endif
            do k=1,levlen(mype)
              do i=1,lonfp
                mn=mn+1
                tmpsnd(mn) = a(i+ioff     ,k+kp,j)
                mn=mn+1
                tmpsnd(mn) = a(i+ioff+lonff,k+kp,j)
              enddo
            enddo
          enddo
        enddo
        lensnd(ii)=mn-locsnd(ii)
      enddo
c
#ifdef REDUCE_GRID
      latlon=0
      do j=1,latlen(mype)
        latlon=latlon+lonfdp(j,mype)
      enddo
#else
      latlon=latlen(mype)*lonlen(mype)
#endif
      mn=0
      do ii=1,ncol
        locrcv(ii)=mn
        ipe=(ii-1)+ipe0
        lenrcv(ii)=levlen(ipe)*latlon*2*nvar
        mn=mn+lenrcv(ii)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   comm_row,ierr)
c restore to k
      mn=0
#ifndef REDUCE_GRID
      lonfp=lonlen(mype)
#endif
      do ii=1,ncol
        ipe=(ii-1)+ipe0
        koff=levstr(ipe)-1
        do n=1,nvar
          kf=(n-1)*levs+levb-1+koff
          do j=1,latlen(mype)
#ifdef REDUCE_GRID
            lonfp=lonfdp(j,mype)
#endif
            do k=1,levlen(ipe)
              do i=1,lonfp
                mn=mn+1
                b(i      ,k+kf,j)=tmprcv(mn)
                mn=mn+1
                b(i+lonfp,k+kf,j)=tmprcv(mn)
              enddo
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
