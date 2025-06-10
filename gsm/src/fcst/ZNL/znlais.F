#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
      subroutine znlais(lons2,a,ist,fac,kd,lat)                                       
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
#include <comznl.h>
      dimension a(lons2)                                                       
      dimension ist(lons2)                                                     
c                                                                               
      nlon=lons2/2                                                             
      do k=1,nst                                                                
        zds(1,k,kd,lat)=0.                                                      
        zds(2,k,kd,lat)=0.                                                      
      enddo                                                                     
      do i=1,nlon                                                               
        i1=i                                                                    
        i2=i+nlon                                                               
        zds(1,1,kd,lat)=zds(1,1,kd,lat)+a(i1)*fac                               
        zds(2,1,kd,lat)=zds(2,1,kd,lat)+a(i2)*fac                               
      enddo                                                                     
      do i=1,nlon                                                               
        i1=i                                                                    
        i2=i+nlon                                                               
        zds(1,ist(i1),kd,lat)=zds(1,ist(i1),kd,lat)+a(i1)*fac                   
        zds(2,ist(i2),kd,lat)=zds(2,ist(i2),kd,lat)+a(i2)*fac                   
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
