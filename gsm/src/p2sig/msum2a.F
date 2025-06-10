#define JCAP jcap
#define JCAP1 jcap1
#define LNT2 lnt2
#define LONF lonf
#define LONF2 lonf2
#define TWOJ1 twoj1
      subroutine msum2a(am,ap,bm,bp,fln,qlnwcs,qder,n)                          
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    msu22       computes negative of divergence operator.          
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-11                     
c                                                                               
c abstract: computes negative of divergence operator for                        
c   the terms a and b from eq 27 (sela 1980, 1982).                             
c   input arrays am, ap, bm, bp consist of combinations                         
c   of northern and southern hemisphere fourier coefficients                    
c   previously computed in symns.                                               
c                                                                               
c program history log:                                                          
c   88-04-11  joseph sela                                                       
c                                                                               
c usage:    call msu22 (am, ap, bm, bp, fln, qlnwcs, qder, n)                   
c   input argument list:                                                        
c     am       - array of multipliers of antisymmetric pln.                     
c     ap       - array of multipliers of     symmetric pln.                     
c     bm       - array of multipliers of antisymmetric pln.                     
c     bp       - array of multipliers of     symmetric pln.                     
c     fln      - array of coefficients which already have                       
c                contributions from other latitudes.                            
c     qlnwcs   - array of legendres times weights.                              
c     qder     - array of derivatives of legendres.                             
c     n        - number of fields in arrays am, ap, bm, bp, fln.                
c                                                                               
c   output argument list:                                                       
c     fln      - array of coefficients modified by the contribution             
c                at the current latitude.                                       
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      dimension am(lonf_,n), ap(lonf2_,n),                                      
     1          bm(lonf_,n), bp(lonf2_,n),                                      
     2 qlnwcs(lnt2_), qder(lnt2_), fln(lnt2_,n)                                 
      dimension       s(lnt2_), sev(twoj1_), sod(twoj1_)                        
cc                                                                              
      npair = (jcap1_-3)/2                                                      
c
      do 2 k=1,n                                                                
      do 220 i=1,twoj1_                                                         
      s(i) = bm(i,k) * qder(i)                                                  
  220 continue                                                                  
      len = twoj1_ - 2                                                          
      do 230 i=1,len                                                            
      s(i+twoj1_) = bp(i,k) * qder(i+twoj1_)                                    
  230 continue                                                                  
      iplus = twoj1_*2 - 2                                                      
      len   = twoj1_ - 4                                                        
cc                                                                              
      do 260 j=1,npair                                                          
      do 240 i=1,len                                                            
      s(i+iplus) = bm(i,k) * qder(i+iplus)                                      
  240 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
cc                                                                              
      do 250 i=1,len                                                            
      s(i+iplus) = bp(i,k) * qder(i+iplus)                                      
  250 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  260 continue                                                                  
cc                                                                              
      do 270 i=1,len                                                            
      s(i+iplus) = bm(i,k) * qder(i+iplus)                                      
  270 continue                                                                  
cc                                                                              
      do 280 i=1,lnt2_                                                          
      fln(i,k) = fln(i,k) + s(i)                                                
  280 continue                                                                  
cc                                                                              
      do 320 l=1,twoj1_,2                                                       
      sev(l  ) =  ap(l+1,k)                                                     
      sev(l+1) = -ap(l  ,k)                                                     
  320 continue                                                                  
      len = twoj1_ - 2                                                          
      do 330 l=1,len,2                                                          
      sod(l  ) =  am(l+1,k)                                                     
      sod(l+1) = -am(l  ,k)                                                     
  330 continue                                                                  
cc                                                                              
      do 420 i=1,twoj1_                                                         
      s(i) = sev(i) * qlnwcs(i)                                                 
  420 continue                                                                  
      len = twoj1_ - 2                                                          
      do 430 i=1,len                                                            
      s(i+twoj1_) = sod(i) * qlnwcs(i+twoj1_)                                   
  430 continue                                                                  
      iplus = twoj1_*2 - 2                                                      
      len   = twoj1_ - 4                                                        
cc                                                                              
      do 460 j=1,npair                                                          
      do 440 i=1,len                                                            
      s(i+iplus) = sev(i) * qlnwcs(i+iplus)                                     
  440 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
cc                                                                              
      do 450 i=1,len                                                            
      s(i+iplus) = sod(i) * qlnwcs(i+iplus)                                     
  450 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  460 continue                                                                  
cc                                                                              
      do 470 i=1,len                                                            
      s(i+iplus) = sev(i) * qlnwcs(i+iplus)                                     
  470 continue                                                                  
cc                                                                              
      do 480 i=1,lnt2_                                                          
      fln(i,k) = fln(i,k) + s(i)                                                
  480 continue                                                                  
    2 continue                                                                  
      return                                                                    
      end                                                                       
