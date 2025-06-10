#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATG2 latg2
#define LATGP latgp
#define LATG2P latg2p

      subroutine preduceg
c
c subroutine preduceg		programmer: hann-ming juang
c
c purpose: partial reduce grid routine for mpi to compute lcapd,lonfd
c          lcapdp and lonfdp for each node.
c
c input:
c	lcapd	(latg2) maximal wave number needed for all latitude
c	lonfd	(latg2) maximal grid number for all latitudes
c	latg2	latitude dimension
c output: all into commpi
c	lonfd	(latg2) full lonf for reduced grid as the same as lonfd
c	lcapdp	(latg2 ,npes) partial wave resolution for reduced grid
c	lonfdp	(latg2p,npes) partial lonf length for reduced grid
c	lonfds	(latg2p,npes) partial lonf start point for reduced grid
c
#include <paramodel.h>
#include <commpi.h>
#ifdef REDUCE_GRID
#include <comreduce.h>
c
      integer, allocatable :: lentmp(:)
      allocate(lentmp(ncol))
c
      latg2=latg2_
      do n=0,npes-1
        llstr=lwvstr(n)
        llens=lwvlen(n)
        do lat=1,latg2
          ls=0
          do ll=1,llens
            l=lwvdef(llstr+ll)
            if( l.lt.lcapd(lat) ) ls=ls+1
          enddo
          lcapdp(lat,n)=ls
        enddo
      enddo
c
      do n=0,nrow-1
        jpe=n*ncol
        jjstr=latstr(jpe)
        jjens=latlen(jpe)
        do lat=1,jjens
          j=latdef(jjstr+lat-1)
          call equdiv(lonfd(j),ncol,lentmp)
          mk=1
          do m=0,ncol-1
            np=n*ncol+m
            lonfds(lat,np)=mk
            lonfdp(lat,np)=lentmp(m+1)
            mk=mk+lonfdp(lat,np)
          enddo
        enddo
      enddo
c
#ifdef DBG
      if(mype.eq.master) then
        do j=1,latg2_
          print *,' lat lcapd lonfd ',j,lcapd(j),lonfd(j)
        enddo
        do n=0,npes-1
          print *,' ----- n= ',n,' ---------'
          do j=1,latg2
            print *,' lat lcapdp ',j,lcapdp(j,n)
          enddo
          do j=1,latlen(n)
            jj=latstr(n)+j-1
            print *,' lat lonfds lonfdp ',jj,lonfds(j,n),lonfdp(j,n)
          enddo
        enddo
      endif
#endif
c
#endif
      return
      end

