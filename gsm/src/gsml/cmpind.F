#define JCAP jcap
      subroutine cmpind                                                         
c                                                                               
#include <paramodel.h>
      parameter(mwave=jcap_)                                                    
c                                                                               
      parameter(mwavep=mwave+1,                                                 
     1          mdim=(mwave+1)*(mwave+2))                                       
c                                                                               
      common/comind/ indxnn(mdim),indxmm(mdim)                                  
c                                                                               
c  indxnn(mdim) :  1-d index of converting input form spher coeff array         
c                  to transposed form array                                     
c  indxmm(mdim) :  1-d index of converting transposed form spher coeff          
c                  array to input form spherical coeff array                    
c                                                                               
      l=0                                                                       
      do 10 m=1,mwavep                                                          
      nend=mwavep-m+1                                                           
      do 10 nn=1,nend                                                           
      n=nn+m-1                                                                  
      l=l+2                                                                     
      indx=(mwavep*(n-m)-(n-m)*(n-m-1)/2+m)*2-1                                 
      indxnn(l-1)=indx                                                          
      indxnn(l  )=indx+1                                                        
   10 continue                                                                  
c                                                                               
      l=0                                                                       
      do 20 nn=1,mwavep                                                         
      lln=mwavep-nn+1                                                           
      do 20 ll=1,lln                                                            
      n=ll+nn-1                                                                 
      m=ll                                                                      
      indx=(m*mwavep-(mwavep-n)-(m-1)*m/2)*2-1                                  
      l=l+2                                                                     
      indxmm(l-1)=indx                                                          
      indxmm(l  )=indx+1                                                        
   20 continue                                                                  
      return                                                                    
c                                                                               
      end                                                                       
