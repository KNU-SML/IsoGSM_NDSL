#define JCAP jcap
#define JCAP1 jcap1
#define JCAP2 jcap2
#define LATG latg
#define LATG2 latg2
#define LNT2 lnt2
#define LNUT lnut
#define LNUT2 lnut2
#define TWOJ1 twoj1
      subroutine pln2t (qlnt, qlnv, colrad, lat)                                
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    pln2t       evaluates associated legendre functions.           
c   prgmmr: joseph sela      org: w/nmc23    date: 88-10-25                     
c                                                                               
c abstract: evaluates the required values of the normalized                     
c   associated legendre function at a prescribed colatitude.                    
c   a standard recursion relation is used with real arithmetic.                 
c                                                                               
c program history log:                                                          
c   88-10-25  joseph sela                                                       
c                                                                               
c usage:    call pln2t (qlnt, qlnv, colrad, lat)                                
c   input argument list:                                                        
c     colrad   - half precision colatitudes in radians for which                
c                the associated legendre functions are to be                    
c                computed.                                                      
c     lat      - index which indicates the current latitude.                    
c                                                                               
c   output argument list:                                                       
c     qlnt     - doubled scalar triangle of                                     
c                half precision associated legendre functions.                  
c     qlnv     - doubled vector triangle of                                     
c                half precision associated legendre functions.                  
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
cra       save                                                                 
cc                                                                              
      dimension         qlnt(lnt2_)                                             
      dimension         qlnv(lnut2_)                                            
      dimension       colrad(latg2_)                                            
cc                                                                              
      dimension         deps(lnut2_)                                            
      real              deps                                                    
      dimension        rdeps(lnut2_)                                            
      real             rdeps                                                    
      dimension            x(jcap1_)                                            
      real                 x                                                    
      dimension           dx(twoj1_)                                            
      real                dx                                                    
      dimension            y(jcap1_)                                            
      real                 y                                                    
      dimension         dpln(lnut2_)                                            
      real              dpln                                                    
      real              colr                                                    
      real            sinlat                                                    
      real              cos2                                                    
      real              prod                                                    
c     common /pln2tc/ifir,jfir,deps,rdeps,dx,y                                  
      common /pln2tc/deps,rdeps,dx,y,jfir                                       
cc                                                                              
      data ifir /0/                                                             
      if  (ifir .eq. 1)  go to 500                                              
           ifir = 1                                                             
      do 200 ll=1,jcap1_                                                        
             rdeps(ll) = 0.0                                                    
  200 continue                                                                  
      lplus = jcap1_                                                            
      len   = jcap1_                                                            
      do 240 inde=2,jcap2_                                                      
      do 220   ll=1,len                                                         
             l = ll - 1                                                         
             n = l + inde - 1                                                   
             rdeps(ll+lplus) = (n*n - l*l) / (4.0 * n*n - 1.0)                  
  220 continue                                                                  
      lplus = lplus + len                                                       
      len = len - 1                                                             
  240 continue                                                                  
      do 260 i=jcap2_,lnut_                                                     
             rdeps(i) = sqrt(rdeps(i))                                          
  260 continue                                                                  
      do 300 i=1,lnut_                                                          
             deps(2*i-1) = rdeps(i)                                             
             deps(2*i  ) = rdeps(i)                                             
  300 continue                                                                  
      ibegin = twoj1_ + 1                                                       
      do 320 i=ibegin,lnut2_                                                    
             rdeps(i) = 1.0/deps(i)                                             
  320 continue                                                                  
      do 400 ll=1,jcap1_                                                        
             x(ll) = ll*2+1                                                     
  400 continue                                                                  
      do 420 ll=1,jcap1_                                                        
             y(ll) = x(ll)/(x(ll)-1.)                                           
  420 continue                                                                  
      do 440 ll=1,jcap1_                                                        
             x(ll) = sqrt(x(ll))                                                
  440 continue                                                                  
      do 460 ll=1,jcap1_                                                        
             dx(2*ll-1) = x(ll)                                                 
             dx(2*ll  ) = x(ll)                                                 
  460 continue                                                                  
  500 continue                                                                  
cc                                                                              
      colr   = colrad(lat)                                                      
      sinlat = cos(colr)                                                        
      cos2   = 1.0 - sinlat * sinlat                                            
      prod   = 1.0                                                              
      do 600 ll=1,jcap1_                                                        
             x(ll) = 0.5*prod                                                   
cccc         if (prod .lt. 1.0e-75)  prod=0.0                                   
             prod = prod*cos2*y(ll)                                             
  600 continue                                                                  
      do 620 ll=1,jcap1_                                                        
             x(ll) = sqrt(x(ll))                                                
  620 continue                                                                  
      do 640 ll=1,jcap1_                                                        
             dpln(2*ll-1) = x(ll)                                               
             dpln(2*ll  ) = x(ll)                                               
  640 continue                                                                  
      lplus = twoj1_                                                            
      do 700 ll=1,twoj1_                                                        
             dpln(ll+lplus) = dx(ll) * sinlat * dpln(ll)                        
  700 continue                                                                  
      lp2 = 0                                                                   
      lp1 =     twoj1_                                                          
      lp0 = 2 * twoj1_                                                          
      len =     twoj1_ - 2                                                      
      do 740  n=3,jcap2_                                                        
cdir$ ivdep                                                                     
      do 720 ll=1,len                                                           
             dpln(ll+lp0) = (sinlat * dpln(ll+lp1)                              
     1               - deps(ll+lp1) * dpln(ll+lp2)) * rdeps(ll+lp0)             
  720 continue                                                                  
      lp2 = lp1                                                                 
      lp1 = lp0                                                                 
      lp0 = lp0 + len                                                           
      len = len - 2                                                             
  740 continue                                                                  
cc                                                                              
      do 800 i=1,lnut2_                                                         
             qlnv(i) = dpln(i)                                                  
  800 continue                                                                  
cc                                                                              
      len = 2 * twoj1_ - 2                                                      
      do 820 ll=1,len                                                           
             qlnt(ll) = qlnv(ll)                                                
  820 continue                                                                  
      lpv = 2 * twoj1_                                                          
      lpt = lpv - 2                                                             
      len = twoj1_ - 4                                                          
      do 860 n=3,jcap1_                                                         
      do 840 ll=1,len                                                           
             qlnt(ll+lpt) = qlnv(ll+lpv)                                        
  840 continue                                                                  
      lpv = lpv + len + 2                                                       
      lpt = lpt + len                                                           
      len = len - 2                                                             
  860 continue                                                                  
cc                                                                              
      return                                                                    
      end                                                                       
