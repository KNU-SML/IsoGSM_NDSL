#define LATG2 latg2
#define LATG latg
#define LATGP latgp
#define LONF lonf
#define LONF2 lonf2
#define LONF2P lonf2p
#define LATG2P latg2p
#define LONFP lonfp
#include <define.h>
      subroutine filtr1(tem,te,dim,di,zem,ze,
     1                  rm,rq,
     1                  qm,q,
     1                  filta,lnt22,lnt2,levs,levh)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
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
     1   q(lnt22),      qm(lnt22),
#ifdef NISLQ
     1  rq(LONF2S,levh,LATG2S), rm(LONF2S,levh,LATG2S) 
#else
     1  rq(lnt22,levh), rm(lnt22,levh)                   
#endif
      filtb = (1.-filta)*0.5                       
      do j=1,lnt2   
        qm(j)=q(j)
      enddo
      do 4900 k=1,levs    
      do 4900 j=1,lnt2   
      tem(j,k)=filtb*tem(j,k)+filta*te(j,k)  
      dim(j,k)=filtb*dim(j,k)+filta*di(j,k)  
      zem(j,k)=filtb*zem(j,k)+filta*ze(j,k)  
4900  continue                               

#ifdef NISLQ
      do j = 1, LATG2S
        do k =1, levh
          do i = 1, LONF2S
            rm(i,k,j)=filtb*rm(i,k,j)+filta*rq(i,k,j)
          enddo
        enddo
      enddo
#else
      do 4901 k=1,levh                                            
      do 4901 j=1,lnt2                                           
        rm(j,k)=filtb*rm(j,k)+filta*rq(j,k)                     
4901  continue                                                     
#endif
      return                                                      
      end                                                        
