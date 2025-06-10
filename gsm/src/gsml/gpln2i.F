#define JCAP1 jcap1
#define JCAP2 jcap2
#define JCAP jcap
#define LNUT2 lnut2
#define LNUT lnut
#define TWOJ1 twoj1
      subroutine gpln2i                                                         
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    gpln2t      sets common for subroutine pln2t.                  
c   prgmmr: joseph sela      org: w/nmc23    date: 91-03-14                     
c                                                                               
c abstract: initializes the constant variables and arrays                       
c   of a common for subroutine pln2t.                                           
c                                                                               
c program history log:                                                          
c   91-03-14  joseph sela                                                       
c                                                                               
c usage:    call gpln2t                                                         
c                                                                               
c remarks: call subroutine once before calls to pln2t.                          
c          refer to pln2t for additional documentation.                         
c                                                                               
c attributes:                                                                   
c   language: fortran, cft77.                                                   
c   machine:  cray y-mp.                                                        
c                                                                               
c$$$                                                                            
cc                                                                              
#include <paramodel.h>
      dimension       x(jcap1_)                                                 
cc                                                                              
      common /pln2ti/ deps(lnut2_),rdeps(lnut2_),                               
     1                dx(twoj1_),y(jcap1_),indxmv(lnut2_)                       
cc                                                                              
ccc   data ifir /0/                                                             
cc                                                                              
ccc   if  (ifir .eq. 1)  go to 500                                              
ccc        ifir = 1                                                             
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
c 500 continue                                                                  
cc                                                                              
cc    set index array for transposing vector array                              
cc    from cray order to ibm order.                                             
      l=0                                                                       
      do 640 nn=1,jcap2_                                                        
      lln=min0(jcap2_-nn+1,jcap1_)                                              
      do 620 ll=1,lln                                                           
      indx=((jcap_+3)*(ll-1)-(ll-1)*ll/2+nn)*2                                  
      l=l+2                                                                     
      indxmv(l-1)=indx-1                                                        
      indxmv(l  )=indx                                                          
  620 continue                                                                  
  640 continue                                                                  
cc                                                                              
      return                                                                    
      end                                                                       
