#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
      subroutine zerznl                                                         
#include <paramodel.h>
#include <comznl.h>
#ifdef MP
#define LATG2S latg2p_
#else
#define LATG2S latg2_
#endif
      do m=1,LATG2S
      do l=1,nrm                                                                
      do k=1,levs_                                                              
      do j=1,2                                                                  
      zdm(j,k,l,m)=0.                                                           
      enddo                                                                     
      enddo                                                                     
      enddo                                                                     
      enddo                                                                     
      do m=1,LATG2S                                                             
      do l=1,nrs                                                                
      do k=1,nst                                                                
      do j=1,2                                                                  
      zds(j,k,l,m)=0.                                                           
      enddo                                                                     
      enddo                                                                     
      enddo                                                                     
      enddo                                                                     
      return                                                                    
      end                                                                       
