#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
      subroutine znlaim(lons2,nlonx,a,ps,kd,lat)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
#include <comznl.h>
      dimension a(nlonx,levs_),ps(nlonx)                                       
c                                                                               
      nlon=lons2/2                                                             
      do 20 k=1,levs_                                                           
      zdm(1,k,kd,lat)=0.                                                        
      zdm(2,k,kd,lat)=0.                                                        
      do 10 i=1,nlon                                                            
      i1=i                                                                      
      i2=i+nlon                                                                 
      zdm(1,k,kd,lat)=zdm(1,k,kd,lat)+a(i1,k)*ps(i1)                            
      zdm(2,k,kd,lat)=zdm(2,k,kd,lat)+a(i2,k)*ps(i2)                            
10    continue                                                                  
20    continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
