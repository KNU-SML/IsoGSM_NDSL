#define JCAP jcap
      subroutine transo(a,kmax)                                                 
c                                                                               
#include <paramodel.h>
      parameter(mwave=jcap_)                                                    
c                                                                               
      parameter(mwavep=mwave+1,                                                 
     1          mdim=(mwave+1)*(mwave+2))                                       
c                                                                               
      common/comind/ indxnn(mdim),indxmm(mdim)                                  
c                                                                               
      dimension a(mdim,kmax)                                                    
      dimension b(mdim)                                                         
c                                                                               
      do 10 k=1,kmax                                                            
      do 11 m=1,mdim                                                            
      b(indxmm(m))=a(m,k)                                                       
   11 continue                                                                  
      do 12 m=1,mdim                                                            
      a(m,k)=b(m)                                                               
   12 continue                                                                  
   10 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
