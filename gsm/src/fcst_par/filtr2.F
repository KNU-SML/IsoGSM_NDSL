#define LATG2 latg2
#define LATG latg
#define LATGP latgp
#define LONF lonf
#define LONFP lonfp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LATG2P latg2p
#include <define.h>
      subroutine filtr2(tem,te,dim,di,zem,ze,rm,rq,
     1                                 filta,lnt22,lnt2,levs,levh)     

#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#ifdef NISLQ
#include <paramodel.h>
#include <constant.h>
#endif
      dimension                                                       
     1  te(lnt22,levs), di(lnt22,levs), ze(lnt22,levs),              
     1 tem(lnt22,levs),dim(lnt22,levs),zem(lnt22,levs),             
#ifdef NISLQ
     1  rq(LONF2S,levh,LATG2S), rm(LONF2S,levh,LATG2S)
#else
     1  rq(lnt22,levh), rm(lnt22,levh)
#endif
      filtb = (1.-filta)* 0.5                                      
      do 4900 k=1,levs                                            
      do 4900 j=1,lnt2                                           
      tem(j,k)=tem(j,k)+filtb*te(j,k)                           
      dim(j,k)=dim(j,k)+filtb*di(j,k)                          
      zem(j,k)=zem(j,k)+filtb*ze(j,k)                         
4900  continue                                               


#ifdef NISLQ
      do j = 1, LATG2S
        do k = 1, levh
          do i = 1, LONF2S
            rm(i,k,j)=rm(i,k,j)+filtb*rq(i,k,j)
          enddo
        enddo
      enddo
#else
      do 4901 k=1,levh                                                
      do 4901 j=1,lnt2                                               
      rm(j,k)=rm(j,k)+filtb*rq(j,k)                                 
4901  continue                                                     
#endif
      return                                                      
      end                                                        
