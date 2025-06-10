#define JCAP jcap
#define JCAP1 jcap1
#define JCAP2 jcap2
#define LNT2 lnt2
#define LNUT2 lnut2
#define RERTH rerth
      subroutine dztouv(di,ze,uu,vv,e)                                          
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    dztouv      computes wind from vorticity, divergence.          
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-13                     
c                                                                               
c abstract: calculates   spectral representation of cosine-weighted             
c   wind components from spectral representation of vorticity and               
c   divergence.                                                                 
c                                                                               
c program history log:                                                          
c   88-05-13  joseph sela                                                       
c                                                                               
c usage:    call dztouv (di, ze, uu, vv, e)                                     
c   input argument list:                                                        
c     di       - triangular spectral divergence at one level.                   
c     ze       - triangular spectral vorticity  at one level.                   
c     e        - array of constants computed in epslon.  used to                
c                calculate cosine-weighted latitudinal derivatives              
c                of stream function or velocity potential.                      
c                                                                               
c   output argument list:                                                       
c     uu       - triangular spectral zonal pseudo-wind at one level.            
c     vv       - triangular spectral meridional pseudo-wind                     
c                at one level.                                                  
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      save                                                                      
      common /scrtch/d(2,jcap1_,jcap1_),z(2,jcap1_,jcap1_),                     
     1   u(2,jcap1_,jcap2_),v(2,jcap1_,jcap2_)                                  
      dimension e(jcap1_,jcap2_)                                                
      dimension di(lnt2_),ze(lnt2_),uu(lnut2_),vv(lnut2_)                       
c       complex cil                                                             
      rjcap1=jcap1_                                                             
c                                                                               
      call shufls(di,d,-1)                                                      
      call shufls(ze,z,-1)                                                      
c       l=0                                                                     
      xn = 0. e0                                                               
c       u(x,1,1) =  e(2,1)*z(y,2,1)                                             
        u(1,1,1) =  e(1,2)*z(1,1,2)                                             
        u(2,1,1) =  e(1,2)*z(2,1,2)                                             
c       v(x,1,1) = -e(2,1)*d(y,2,1)                                             
        v(1,1,1) = -e(1,2)*d(1,1,2)                                             
        v(2,1,1) = -e(1,2)*d(2,1,2)                                             
      do 1 i=2,jcap_                                                            
      xn = xn + 1. e0                                                          
      r = 1. e0/xn                                                             
      r1 = 1. e0/(xn+1. e0)                                                   
c     u(x,i,1)=-e(i,1)*z(y,i-1,1)*r+e(i+1,1)*z(y,i+1,1)*r1                      
      u(1,1,i)=-e(1,i)*z(1,1,i-1)*r+e(1,i+1)*z(1,1,i+1)*r1                      
      u(2,1,i)=-e(1,i)*z(2,1,i-1)*r+e(1,i+1)*z(2,1,i+1)*r1                      
c     v(x,i,1)= e(i,1)*d(y,i-1,1)*r-e(i+1,1)*d(y,i+1,1)*r1                      
      v(1,1,i)= e(1,i)*d(1,1,i-1)*r-e(1,i+1)*d(1,1,i+1)*r1                      
      v(2,1,i)= e(1,i)*d(2,1,i-1)*r-e(1,i+1)*d(2,1,i+1)*r1                      
    1 continue                                                                  
      xn = xn + 1. e0                                                          
      r = 1. e0/xn                                                             
c       u(x,jcap1_,1) = -e(jcap1_,1)*z(y,jcap_,1)*r                             
        u(1,1,jcap1_) = -e(1,jcap1_)*z(1,1,jcap_)*r                             
        u(2,1,jcap1_) = -e(1,jcap1_)*z(2,1,jcap_)*r                             
c       v(x,jcap1_,1) =  e(jcap1_,1)*d(y,jcap_,1)*r                             
        v(1,1,jcap1_) =  e(1,jcap1_)*d(1,1,jcap_)*r                             
        v(2,1,jcap1_) =  e(1,jcap1_)*d(2,1,jcap_)*r                             
      xn = xn + 1. e0                                                          
      r = 1. e0/xn                                                             
c       u(x,jcap2_,1) = -e(jcap2_,1)*z(y,jcap1_)*r                              
        u(1,1,jcap2_) = -e(1,jcap2_)*z(1,1,jcap1_)*r                            
        u(2,1,jcap2_) = -e(1,jcap2_)*z(2,1,jcap1_)*r                            
c       v(x,jcap2_,1) =  e(jcap2_,1)*d(y,jcap1_)*r                              
        v(1,1,jcap2_) =  e(1,jcap2_)*d(1,1,jcap1_)*r                            
        v(2,1,jcap2_) =  e(1,jcap2_)*d(2,1,jcap1_)*r                            
c         fin l=0                                                               
      do 1000 l=2,jcap1_                                                        
      xl = l                                                                    
      xn = xl-1. e0                                                            
      xll = xn                                                                  
      r = 1. e0/xn                                                             
      r1 = 1. e0/xl                                                            
c       cil=cmplx(0.,1.)*(l-1)                                                  
c       u(x,1,l)=-cil*d(y,1,l)*r*r1 + e(2,l)*z(y,2,l)*r1                        
        u(1,l,1)= d(2,l,1)*r*r1*xll+e(l,2)*z(1,l,2)*r1                          
        u(2,l,1)=-d(1,l,1)*r*r1*xll+e(l,2)*z(2,l,2)*r1                          
c       v(x,1,l)=-cil*z(y,1,l)*r*r1 - e(2,l)*d(y,2,l)*r1                        
        v(1,l,1)= z(2,l,1)*r*r1*xll-e(l,2)*d(1,l,2)*r1                          
        v(2,l,1)=-z(1,l,1)*r*r1*xll-e(l,2)*d(2,l,2)*r1                          
      do 2 i=2,jcap_                                                            
      xn = xn+1. e0                                                            
      if(xn.gt.rjcap1) go to 2                                                  
      r = 1. e0/xn                                                             
      r1 = 1. e0/(xn+1. e0)                                                   
      r2 = r*r1                                                                 
      rl = r2*xll                                                               
c       u(x,i,l)=-e(i,l)*z(y,i-1,l)*r-cil*d(y,i,l)*r2                           
c                      + e(i+1,l)*z(y,i+1,l)*r1                                 
      u(1,l,i)=-e(l,i)*z(1,l,i-1)*r+d(2,l,i)*rl+e(l,i+1)*z(1,l,i+1)*r1          
      u(2,l,i)=-e(l,i)*z(2,l,i-1)*r-d(1,l,i)*rl+e(l,i+1)*z(2,l,i+1)*r1          
c       v(x,i,l)= e(i,l)*d(y,i-1,l)*r-cil*z(y,i,l)*r2                           
c                      - e(i+1,l)*d(y,i+1,l)*r1                                 
      v(1,l,i)= e(l,i)*d(1,l,i-1)*r+z(2,l,i)*rl-e(l,i+1)*d(1,l,i+1)*r1          
      v(2,l,i)= e(l,i)*d(2,l,i-1)*r-z(1,l,i)*rl-e(l,i+1)*d(2,l,i+1)*r1          
    2 continue                                                                  
      xn = xn +1. e0                                                           
      if(xn.gt.rjcap1) go to 1000                                               
      r = 1. e0/xn                                                             
      r1 = r/(xn+1. e0)                                                        
c       u(x,jcap1_,l)=-e(jcap1_,l)*z(y,jcap_,l)*r-cil*d(y,jcap1_,l)*r1          
      u(1,l,jcap1_)=-e(l,jcap1_)*z(1,l,jcap_)*r+d(2,l,jcap1_)*r1*xll            
      u(2,l,jcap1_)=-e(l,jcap1_)*z(2,l,jcap_)*r-d(1,l,jcap1_)*r1*xll            
c       v(x,jcap1_,l)= e(jcap1_,l)*d(y,jcap_,l)*r-cil*z(y,jcap1_,l)*r1          
      v(1,l,jcap1_)= e(l,jcap1_)*d(1,l,jcap_)*r+z(2,l,jcap1_)*r1*xll            
      v(2,l,jcap1_)= e(l,jcap1_)*d(2,l,jcap_)*r-z(1,l,jcap1_)*r1*xll            
      xn = xn + 1. e0                                                          
      r = 1. e0/xn                                                             
c       u(x,jcap2_,l) = -e(jcap2_,l)*z(y,jcap1_,l)*r                            
        u(1,l,jcap2_) = -e(l,jcap2_)*z(1,l,jcap1_)*r                            
        u(2,l,jcap2_) = -e(l,jcap2_)*z(2,l,jcap1_)*r                            
c       v(x,jcap2_,l) =  e(jcap2_,l)*d(y,jcap1_,l)*r                            
        v(1,l,jcap2_) =  e(l,jcap2_)*d(1,l,jcap1_)*r                            
        v(2,l,jcap2_) =  e(l,jcap2_)*d(2,l,jcap1_)*r                            
 1000 continue                                                                  
c....                                                                           
      call shuflv(uu,u,1)                                                       
      call shuflv(vv,v,1)                                                       
      do 6 j=1,lnut2_                                                           
      uu(j)=uu(j)*rerth_                                                        
      vv(j)=vv(j)*rerth_                                                        
    6 continue                                                                  
      return                                                                    
      end                                                                       
