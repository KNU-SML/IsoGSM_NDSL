#define JCAP jcap
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF lonf
#define NPES npes
      subroutine mpgpfk2fpk(a,lonf2p,latg2p,ntotal,
     1                      b,lonf2 ,latg2 ,nlenmx,
     2                      nstr,nend)
c
c  mpi communication from a(ip,jp,kf) to b(if,jf,kp) using npes tasks
c
      implicit none
#include <paramodel.h>
#include <commpi.h>
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
      integer lonf2p,latg2p,lonf2,latg2,ntotal,nlenmx
      integer ijk,n,k,j,i,ii,jj,kk,ierr
      integer lonlend,lonstrd,lonff
      real    a,b
      dimension a(lonf2p,latg2p,ntotal),b(lonf2,latg2,nlenmx)
      integer nstr(0:npes-1),nend(0:npes-1)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::lensnd(:),locsnd(:),lenrcv(:),locrcv(:)
c
      allocate(tmpsnd(lonf2p*latg2p*ntotal))
      allocate(tmprcv(lonf2 *latg2 *nlenmx))
      allocate(lensnd(0:npes_-1),lenrcv(0:npes_-1))
      allocate(locsnd(0:npes_-1),locrcv(0:npes_-1))
c
c  slice (ip,jp,kf) in vertical into npes pieces
c
      ijk=0
      do n=0,npes-1
        locsnd(n)=ijk
        do k=nstr(n),nend(n)
          do j=1,latlen(mype)
#ifdef REDUCE_GRID
            lonlend=lonfdp(j,mype)
#else
            lonlend=lonlen(mype)
#endif
            do i=1,lonlend
              ijk=ijk+1
              tmpsnd(ijk)=a(i,j,k)
              ijk=ijk+1
              tmpsnd(ijk)=a(i+lonlend,j,k)
            enddo
          enddo
        enddo
        lensnd(n)=ijk-locsnd(n)
      enddo
c
c  prepare location to receive from all others
c
      ijk=0
      do n=0,npes-1
        locrcv(n)=ijk
#ifdef REDUCE_GRID
        jj=0
        do j=1,latlen(n)
          jj=jj+lonfdp(j,n)
        enddo
        lenrcv(n)=jj*2*(nend(mype)-nstr(mype)+1)
#else
        lenrcv(n)=lonlen(n)*2*latlen(n)*(nend(mype)-nstr(mype)+1)
#endif
        ijk=ijk+lenrcv(n)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   MPI_COMM_GSM,ierr)
c
      ijk=0
      do n=0,npes-1
#ifndef REDUCE_GRID
        lonlend=lonlen(n)
        lonstrd=lonstr(n)-1
        lonff=lonf_
#endif
        do k=nstr(mype),nend(mype)
          kk=k-nstr(mype)+1
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
              ijk=ijk+1
              b(ii,jj,kk) = tmprcv(ijk)
              ii=ii+lonff
              ijk=ijk+1
              b(ii,jj,kk) = tmprcv(ijk)
            enddo
          enddo
         enddo
      enddo
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(lensnd,lenrcv)
      deallocate(locsnd,locrcv)
c
      return
      end
