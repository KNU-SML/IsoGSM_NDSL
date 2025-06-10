#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
#define LATG2 latg2
#define LATG latg
#define LNT2 lnt2
#define LNUT2 lnut2
#define TWOJ1 twoj1
      subroutine pln2i(qlnt,qlnv,colrad,lat)                                    
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
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
cc                                                                              
#include <paramodel.h>
      dimension         qlnt(lnt2_)                                             
      dimension         qlnv(lnut2_)                                            
      dimension       colrad(latg2_)                                            
cc                                                                              
      dimension            x(jcap1_)                                            
      dimension         dpln(lnut2_)                                            
cc                                                                              
      common /pln2ti/ deps(lnut2_),rdeps(lnut2_),                               
     1                dx(twoj1_),y(jcap1_),indxmv(lnut2_)                       
cc                                                                              
ccc   data ifir /0/                                                             
ccc         part between guards made into sr gpln2i.                            
ccc         7 dec 1990      m. rozwodoski                                       
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
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
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
cc    transpose vector dpln array from cray order to ibm order.                 
      do 800 i=1,lnut2_                                                         
             qlnv(indxmv(i)) = dpln(i)                                          
  800 continue                                                                  
cc                                                                              
      lpv = 0                                                                   
      lpt = 0                                                                   
      len = twoj1_                                                              
      do 860  n=1,jcap1_                                                        
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
