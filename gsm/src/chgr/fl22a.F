#define JCAP jcap
#define JCAP1 jcap1
#define LNT2 lnt2
#define LONF lonf
#define LONF2 lonf2
#define TWOJ1 twoj1
      subroutine fl22a(fp,fm,fln,qln,n)                                         
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    fl22        accumulates spectral from fourier input.           
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-04                     
c                                                                               
c abstract: accumulates spectral coefficients from fourier input.               
c   call to subroutine accumulates the contribution at the given                
c   latitude only.                                                              
c                                                                               
c program history log:                                                          
c   88-04-04  joseph sela                                                       
c                                                                               
c usage:    call fl22 (fp, fm, fln, qln, n)                                     
c   input argument list:                                                        
c     fp       - array of fourier  coefficients for  symmetric part.            
c     fm       - array of fourier  coefficients for asymmetric part.            
c     fln      - array of spectral coefficients to be summed into               
c                at the current latitude.                                       
c     qln      - array of associated legendre functions.                        
c     n        - number of fields in fp, fm, fln arrays.                        
c                                                                               
c   output argument list:                                                       
c     fln      - array of spectral coefficients with contribution               
c                from the current latitude.                                     
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      dimension fp(lonf2_,n), fm(lonf_,n),                                      
     1 qln(lnt2_), fln(lnt2_,n)                                                 
      dimension       s(lnt2_)                                                  
cc                                                                              
      npair = (jcap1_-3)/2                                                      
c
      do 2 k=1,n                                                                
      do 220 i=1,twoj1_                                                         
      s(i) = fp(i,k) * qln(i)                                                   
  220 continue                                                                  
      len = twoj1_ - 2                                                          
      do 230 i=1,len                                                            
      s(i+twoj1_) = fm(i,k) * qln(i+twoj1_)                                     
  230 continue                                                                  
      iplus = twoj1_*2 - 2                                                      
      len   = twoj1_ - 4                                                        
cc                                                                              
      do 260 j=1,npair                                                          
      do 240 i=1,len                                                            
      s(i+iplus) = fp(i,k) * qln(i+iplus)                                       
  240 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
cc                                                                              
      do 250 i=1,len                                                            
      s(i+iplus) = fm(i,k) * qln(i+iplus)                                       
  250 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  260 continue                                                                  
cc                                                                              
      do 270 i=1,len                                                            
      s(i+iplus) = fp(i,k) * qln(i+iplus)                                       
  270 continue                                                                  
cc                                                                              
      do 280 i=1,lnt2_                                                          
      fln(i,k) = fln(i,k) + s(i)                                                
  280 continue                                                                  
    2 continue                                                                  
      return                                                                    
      end                                                                       
