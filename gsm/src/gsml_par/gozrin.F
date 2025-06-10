#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
#define LNT2 lnt2
#define LNUT2 lnut2
#define TWOJ1 twoj1
      subroutine gozrin(qlnt,qlnv,qdert,epsi)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
cc                                                                              
#include <paramodel.h>
      dimension          qlnt(lnt2_)                                            
      dimension          qlnv(lnut2_)                                           
      dimension         qdert(lnt2_)                                            
      dimension          epsi(jcap2_,jcap1_)                                    
cc                                                                              
      common /gozcom/ dxa(lnt2_),dxb(lnt2_)                                     
cc                                                                              
ccc         part between guards made into sr ggozri.                            
ccc         7 dec 1990      m. rozwodoski                                       
cc                                                                              
cc    compute pln derivatives in ibm order.                                     
cc                                                                              
      lp0 = 0                                                                   
      lp1 = 2                                                                   
      len = twoj1_                                                              
      do 640  i=1,jcap1_                                                        
      do 620 ll=1,len                                                           
             qdert(ll+lp0) = qlnv(ll+lp1) * dxb(ll+lp0)                         
  620 continue                                                                  
      lp1 = lp1 + len + 2                                                       
      lp0 = lp0 + len                                                           
      len = len - 2                                                             
  640 continue                                                                  
cc                                                                              
                  lend = lnt2_ - 4                                              
      do 720 ll=1,lend                                                          
             qdert(ll+2) = qdert(ll+2) + qlnt(ll) * dxa(ll+2)                   
  720 continue                                                                  
cc                                                                              
      return                                                                    
      end                                                                       
