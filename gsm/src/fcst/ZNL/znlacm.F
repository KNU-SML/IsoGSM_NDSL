#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine znlacm(lons2,nlonx,a,ps,dt,kd,lat)      
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
#include <comznl.h>
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
      dimension a(nlonx,levs_),ps(LONF2S)
c                                                                               
      nlon=lons2/2                                                             
      do 20 k=1,levs_                                                           
      do 10 i=1,nlon                                                            
      i1=i                                                                      
      i2=i+nlon                                                                 
      zdm(1,k,kd,lat)=zdm(1,k,kd,lat)+a(i1,k)*ps(i1)*dt                         
      zdm(2,k,kd,lat)=zdm(2,k,kd,lat)+a(i2,k)*ps(i2)*dt                         
10    continue                                                                  
20    continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
