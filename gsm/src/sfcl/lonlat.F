#define PI pi
      subroutine lonlat(xlon,xlat,colrad,lon,lat)
#include <paramodel.h>
      parameter (tpi=2.e0*pi_,hpi=0.5e0*pi_)
      dimension xlon(lon,lat),xlat(lon,lat),colrad(lat/2)
c....  get normal lon,lat in radians
      bphi = tpi / lon
      do j=1,lat
        do i=1,lon
          xlon(i,j) = (i-1) * bphi
        enddo
      enddo
      do j=1,lat/2
        do i=1,lon
          xlat(i,j) = hpi - colrad(j)
        enddo
      enddo
      do j=lat/2+1,lat
        jj=lat+1-j
        do i=1,lon
          xlat(i,j) = - xlat(i,jj)
        enddo
      enddo
      return
      end
