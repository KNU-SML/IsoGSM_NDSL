#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
#define LNT22 lnt22
      subroutine bartri (f, g, fgbar)                                           
#include <paramodel.h>
      dimension f(lnt22_),g(lnt22_)                                             
c                                                                               
      joff(n,l)=(jcap1_)*(jcap2_)-(jcap1_-l)*(jcap2_-l)+2*(n-l)                 
c                                                                               
      l=0                                                                       
      fgbar = 0.                                                                
      do 1 n=0, jcap_                                                           
      fgbar = fgbar + f(joff(n,l)+1)*g(joff(n,l)+1)                             
1     continue                                                                  
      do 11 n=0, jcap_                                                          
      fgbar = fgbar + f(joff(n,l)+2)*g(joff(n,l)+2)                             
11    continue                                                                  
      fgbar=fgbar*0.5                                                           
      do 3 l=1, jcap_                                                           
      do 2 n=l, jcap_                                                           
      fgbar = fgbar + f(joff(n,l)+1)*g(joff(n,l)+1)                             
2     continue                                                                  
      do 22 n=l, jcap_                                                          
      fgbar = fgbar + f(joff(n,l)+2)*g(joff(n,l)+2)                             
22    continue                                                                  
3     continue                                                                  
      fgbar = sqrt(fgbar)                                                       
      return                                                                    
      end                                                                       
