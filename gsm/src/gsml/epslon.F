#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
      subroutine epslon(eps,jcap)                                               
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    epslon      computes eps, a function of wave number.           
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-01                     
c                                                                               
c abstract: computes eps, a function of wave number.                            
c   eps is used in calculating legendre polys. and their derivatives.           
c   eps is also used in computing winds from divergence and vorticity.          
c                                                                               
c program history log:                                                          
c   88-04-01  joseph sela                                                       
c                                                                               
c usage:    call epslon (eps, jcap)                                             
c   input argument list:                                                        
c     jcap     - index indicating the spectral truncation used.                 
c                                                                               
c   output argument list:                                                       
c     eps      - array computed from sqrt((n**2-l**2)/(4*n**2-1)).              
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension eps(jcap1_,jcap2_)                                              
      jcap1 = jcap + 1                                                          
      jcap2 = jcap + 2                                                          
      do 1 ll=1,jcap1                                                           
      l = ll - 1                                                                
      do 1 inde=2,jcap2                                                         
      n = l + inde - 1                                                          
      a = (n*n - l*l) / (4.0 * n*n - 1.0)                                       
      eps(ll,inde)=sqrt(a)                                                    
1     continue                                                                  
      do 2 ll=1,jcap1_                                                          
      eps(ll,1) = 0.0e0                                                        
2     continue                                                                  
      return                                                                    
      end                                                                       
