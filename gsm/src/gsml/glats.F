#define LATG2 latg2
#define LATG latg
#define PI pi
      subroutine glats(lgghaf,colrad,wgt,wgtcs,rcs2)                            
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    glats       computes location of gaussian latitudes.           
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-05                     
c                                                                               
c abstract: computes the location of the gaussian latitudes for the             
c   input lgghaf.  the latitudes are determined by finding                      
c   the zeros of the legendre polynomials.                                      
c                                                                               
c program history log:                                                          
c   88-04-05  joseph sela                                                       
c                                                                               
c usage:    call glats (lgghaf, colrad, wgt, wgtcs, rcs2)                       
c   input argument list:                                                        
c     lgghaf   - number of gaussian latitudes in a hemisphere.                  
c                                                                               
c   output argument list:                                                       
c     colrad   - array of colatitude of gaussian latitudes                      
c                in northern hemisphere.                                        
c     wgt      - array of weights at each gaussian latitude                     
c                required for gaussian quadrature.                              
c     wgtcs    - array of gaussian weight/sin of colatitude squared.            
c     rcs2     - array of reciprocal  of  sin of colatitude squared.            
c                                                                               
c   output files:                                                               
c     output   - printout file.                                                 
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
cccc  half precision colrad,wgt,wgtcs,rcs2                                      
#include <paramodel.h>
      real colrad,wgt,wgtcs,rcs2                                              
      dimension colrad( latg2_ ),wgt( latg2_ ),wgtcs( latg2_ )                  
      dimension rcs2( latg2_ )                                                  
      eps=1.e-12                                                                
      si = 1.0                                                                  
      l2=2*lgghaf                                                               
      rl2=l2                                                                    
      scale = 2.0/(rl2*rl2)                                                     
      k1=l2-1                                                                   
      dradz = pi_ / 360.
      rad = 0.0                                                                 
      do 1000 k=1,lgghaf                                                        
      iter=0                                                                    
      drad=dradz                                                                
1     call poly(l2,rad,p2)                                                      
2     p1 =p2                                                                    
      iter=iter+1                                                               
      rad=rad+drad                                                              
      call poly(l2,rad,p2)                                                      
      if(sign(si,p1).eq.sign(si,p2)) go to 2                                    
      if(drad.lt.eps)go to 3                                                    
      rad=rad-drad                                                              
      drad = drad * 0.25                                                        
      go to 1                                                                   
3     continue                                                                  
      colrad(k)=rad                                                             
      phi = rad * 180 / pi_                                                     
      call poly(k1,rad,p1)                                                      
      x = cos(rad)                                                              
      w = scale * (1.0 - x*x)/ (p1*p1)                                          
      wgt(k) = w                                                                
      sn = sin(rad)                                                             
      w=w/(sn*sn)                                                               
      wgtcs(k) = w                                                              
      rc=1./(sn*sn)                                                             
      rcs2(k) = rc                                                              
      call poly(l2,rad,p1)                                                      
1000  continue                                                                  
      return                                                                    
      end                                                                       
