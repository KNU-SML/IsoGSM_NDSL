#define JCAP jcap
#define JCAP2 jcap2
#define LNUT2 lnut2
#define LONF lonf
#define LONF2 lonf2
#define TWOJ1 twoj1
      subroutine sumv2a (vln,ap,qln,n)                                          
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    sumvns      produces vector fourier coefficients.              
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-04                     
c                                                                               
c abstract: uses spherical harmonics of a pseudo vector to produce              
c   the fourier coefficients of the fields on the northern                      
c   and southern hemisphere latitude circles for which                          
c   the input legendre functions apply.                                         
c                                                                               
c program history log:                                                          
c   88-04-04  joseph sela                                                       
c                                                                               
c usage:    call sumvns (vln, ap, qln, n)                                       
c   input argument list:                                                        
c     vln      - array of spherical harmonics.                                  
c     qln      - array of associated legendre functions.                        
c     n        - number of fields in vln and ap arrays.                         
c                                                                               
c   output argument list:                                                       
c     ap       - array of fourier representations of fields on                  
c                northern and southern hemisphere latitude circles.             
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      dimension ap(lonf2_,n),                                                   
     1 qln(lnut2_), vln(lnut2_,n)                                               
      dimension       s(lnut2_), sev(twoj1_), sod(twoj1_)                       
cc                                                                              
      npair = (jcap2_-2)/2                                                      
#ifdef CRAY_THREAD
cmic$ do all private(s,sev,sod) autoscope                                       
#endif
      do 2 k=1,n                                                                
      do 100 i=1,lnut2_                                                         
      s(i) = qln(i) * vln(i,k)                                                  
  100 continue                                                                  
      do 200 i=1,twoj1_                                                         
      sev(i) = s(i)                                                             
      sod(i) = s(i+twoj1_)                                                      
  200 continue                                                                  
      iplus = twoj1_*2                                                          
      len   = twoj1_-2                                                          
cc                                                                              
      do 600 j=1,npair                                                          
      do 400 i=1,len                                                            
      sev(i) = sev(i) + s(i+iplus)                                              
  400 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
cc                                                                              
      do 500 i=1,len                                                            
      sod(i) = sod(i) + s(i+iplus)                                              
  500 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  600 continue                                                                  
cc                                                                              
      do 800 i=1,twoj1_                                                         
      ap(i      ,k) = sev(i) + sod(i)                                           
  800 continue                                                                  
      do 900 i=1,twoj1_                                                         
      ap(i+lonf_,k) = sev(i) - sod(i)                                           
  900 continue                                                                  
    2 continue                                                                  
      return                                                                    
      end                                                                       
