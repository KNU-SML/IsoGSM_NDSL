#include <define.h>
#ifdef SCM
      subroutine sincoslat(sinlat,coslat,xlon,xlat,lon,lat)
#else
      subroutine sincoslat(sinlat,coslat,colrad,lon,lat)
#endif
      dimension sinlat(lon,lat),coslat(lon,lat),colrad(lat/2)                    
#ifdef SCM
      dimension xlat(lon,lat),xlon(lon,lat)
      do j = 1, lat
      do i = 1, lon
        sinlat(i,j) = sin(xlat(i,j))
        coslat(i,j) = sqrt(1. e 0 - sinlat(i,j)*sinlat(i,j))
      enddo
      enddo
#else
c....  get normal sinlat and coslat 
      do j=1,lat/2
        sinlaj = cos(colrad(j))
        coslaj = sqrt(1. e 0 - sinlaj*sinlaj)
        do i=1,lon
          sinlat(i,j) = sinlaj
          coslat(i,j) = coslaj
        enddo
      enddo
      do j=lat/2+1,lat
        jj=lat+1-j
        do i=1,lon
          sinlat(i,j) = -sinlat(i,jj)
          coslat(i,j) =  coslat(i,jj)
        enddo
      enddo
#endif
      return                                                                    
      end                                                                       
