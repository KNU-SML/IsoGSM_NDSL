#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
#define LNT2 lnt2
#define TWOJ1 twoj1
      subroutine ggozrm(epsi)                                                   
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    ggozrm      sets common for subroutine gozrmt.                 
c   prgmmr: joseph sela      org: w/nmc23    date: 91-03-14                     
c                                                                               
c abstract: initializes the constant variables and arrays                       
c   of a common for subroutine gozrmt.                                          
c                                                                               
c program history log:                                                          
c   91-03-14  joseph sela                                                       
c                                                                               
c usage:    call ggozrm ( eps )                                                 
c                                                                               
c remarks: call subroutine once before calls to gozrmt.                         
c          refer to gozrmt for additional documentation.                        
c                                                                               
c attributes:                                                                   
c   language: fortran, cft77.                                                   
c   machine:  cray y-mp.                                                        
c                                                                               
c$$$                                                                            
cc                                                                              
#include <paramodel.h>
      dimension          epsi(jcap2_,jcap1_)                                    
cc                                                                              
      dimension         dxint(lnt2_)                                            
      dimension            dx(twoj1_,jcap2_)                                    
      dimension          deps(twoj1_,jcap2_)                                    
cc                                                                              
      common /gozcom/ dxa(lnt2_),dxb(lnt2_)                                     
cc                                                                              
      do 200 ll=1,twoj1_                                                        
             dxint(2*ll-1) = ll                                                 
             dxint(2*ll  ) = ll                                                 
  200 continue                                                                  
      lp = 0                                                                    
      do 240  i=1,jcap2_                                                        
      do 220 ll=1,twoj1_                                                        
             dx(ll,i) = dxint(ll+lp)                                            
  220 continue                                                                  
      lp = lp + 2                                                               
  240 continue                                                                  
      do 280  i=1,jcap2_                                                        
      do 260 ll=1,jcap1_                                                        
             deps(2*ll-1,i) = epsi(i,ll)                                        
             deps(2*ll  ,i) = epsi(i,ll)                                        
  260 continue                                                                  
  280 continue                                                                  
      do 300 ll=1,twoj1_                                                        
             dxa(ll) = 0.0                                                      
  300 continue                                                                  
      lp1 = twoj1_                                                              
      len = twoj1_ - 2                                                          
      do 340  i=1,jcap_                                                         
      do 320 ll=1,len                                                           
             dxa(ll+lp1) =  dx(ll,i+1) * deps(ll,i+1)                           
             dxb(ll+lp1) = -dx(ll,i  ) * deps(ll,i+2)                           
  320 continue                                                                  
      lp1 = lp1 + len                                                           
      len = len - 2                                                             
  340 continue                                                                  
      do 380  i=1,jcap2_                                                        
      do 360 ll=1,twoj1_                                                        
             dx(ll,i) = dx(ll,i) - 1.e0                                        
  360 continue                                                                  
  380 continue                                                                  
      do 400 ll=1,twoj1_                                                        
             dxb(ll) = -dx(ll,1) * deps(ll,2)                                   
  400 continue                                                                  
cc                                                                              
cc    transpose scalar arrays dxa, dxb, dxc, dxd                                
cc    from cray order to ibm order.                                             
      call transo (dxa, 2)                                                      
cc                                                                              
      return                                                                    
      end                                                                       
