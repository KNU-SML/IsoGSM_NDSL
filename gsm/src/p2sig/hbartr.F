#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine hbartr(f, g, fgbar)                                            
#include <paramodel.h>
      save                                                                      
c                                                                               
#include <paramter.h>
c                                                                               
      parameter(mwavep=mwave+1,                                                 
     1          mdim=(mwave+1)*(mwave+1)*2)                                     
c                                                                               
      dimension f(mdim)                                                         
      dimension g(mdim)                                                         
c                                                                               
      fgbar = 0.                                                                
      l=0                                                                       
      do 1 i=1,mwavep                                                           
      l=l+2                                                                     
      fgbar = fgbar + f(l-1)*g(l-1)*0.5                                         
      nnmax=mwavep-i+1                                                          
      do 2 j=2,nnmax                                                            
      l=l+2                                                                     
      fgbar = fgbar+f(l-1)*g(l-1)+f(l  )*g(l  )                                 
    2 continue                                                                  
    1 continue                                                                  
      fgbar =  sqrt  (fgbar)                                                    
c                                                                               
      return                                                                    
      end                                                                       
