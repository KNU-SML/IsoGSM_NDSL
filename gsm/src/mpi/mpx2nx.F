#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATGP latgp
#define LATG2P latg2p
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
      subroutine mpx2nx(a,lonf2,lota,b,lonf2p,lotb,latg2p,
     &                  leva,levb,nlev)
c$$$  subprogram documentation block
c
c subprogram:   mpx2nx
c
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:   transpose (if,k,jp) to (ip,k,jp)
c
c program history log:
c    99-06-27  henry juang      finish entire test for gsm
c
c usage:        call mpx2nx(a,lonf2,lota,b,lonf2p,lotb,latg2p,
c        &                  leva,levb,nlev)
c
c    input argument lists:
c       a       - real (lonf2,lota,latg2p) partial field in j k
c       lonf2p  - integer partial longitude grid x 2
c       latg2p  - integer partial latitude grid / 2
c       lonf2   - integer total longitude grid x 2
c       leva    - integer starting vertical layer for a
c       levb    - integer starting vertical layer for b
c       nlev    - integer total set of fields
c
c    output argument list:
c       b       - real (lonf2p,lotb,latg2p) partial field in i j
c
c subprograms called:
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
      integer lonf2,lonf2p,latg2p,lota,lotb,leva,levb,nlev
     1       ,i,j,n,ioff,noffa,noffb,lonf2d,lonfp,lonff
      real    a,b
      dimension a(lonf2,lota,latg2p),b(lonf2p,lotb,latg2p)
c
c option for 1-d decomposition
      if( ncol.eq.1 ) then
        do j=1,latg2p
          do n=1,lota
            do i=1,lonf2
              b(i,n,j) = a(i,n,j)
            enddo
          enddo
        enddo
        return
      endif
c
c otherwise
c cut the need and through away the remain.
#ifndef REDUCE_GRID
      lonff=lonf_
      lonfp=lonlen(mype)
      ioff=lonstr(mype)-1
#endif
      do n=1,nlev
        noffa=leva-1+n
        noffb=levb-1+n
        do j=1,latlen(mype)
#ifdef REDUCE_GRID
          lonff=lonfd(latdef(latstr(mype)+j-1))
          lonfp=lonfdp(j,mype)
          ioff=lonfds(j,mype)-1
#endif
          do i=1,lonfp
            b(i      ,noffb,j)=a(i+ioff      ,noffa,j)
            b(i+lonfp,noffb,j)=a(i+ioff+lonff,noffa,j)
          enddo
        enddo
      enddo
c
      return
      end
