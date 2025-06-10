#define PI pi
#include <define.h>
#ifdef RSM
#endif
      subroutine getlatlon(idim,jdim,rlat,rlon)
c
      implicit none
#include <paramodel.h>
      integer idim,jdim
      real rlat(idim,jdim),rlon(idim,jdim)
      real rdelx,rdely,dlamda0,dx
      real gaul(jdim)
c
      integer i,j
c
c  compute latitude and longitude
c
#ifdef RSM
      call setgrd(rlat,rlon,rdelx,rdely,dlamda0)
      call shalfo(rlat,1)
      call shalfo(rlon,1)
c
      do j=1,jdim
        do i=1,idim
          rlat(i,j) = rlat(i,j) * 180. / pi_
          rlon(i,j) = rlon(i,j) * 180. / pi_
        enddo
      enddo
#else
c
c  compute gaussian latitude and longitude
c
      dx=360./float(idim)
      call gaulat(gaul,jdim)
      do j=1,jdim
        do i=1,idim
          rlat(i,j)=90.-gaul(j)
          rlon(i,j)=float(i-1)*dx
          if(rlon(i,j).gt.180.) then
            rlon(i,j)=rlon(i,j)-360.
          endif
        enddo
      enddo
#endif
c
      return
      end
