#define JCAP jcap
#define JCAP1 jcap1
#define JCAP2 jcap2
#define LATG latg
#define LATG2 latg2
#define LNT2 lnt2
#define LNUT2 lnut2
#define RERTH rerth
#define TWOJ1 twoj1
      subroutine gozrmt (qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgtl)              
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    gozrmt      computes derivatives of legendres.                 
c   prgmmr: joseph sela      org: w/nmc23    date: 88-11-02                     
c                                                                               
c abstract: computes derivatives of associated legendre functions               
c   and for convenience, other required products of                             
c   legendres and factors involving wave number and latitude.                   
c   the resulting arrays are required for the application                       
c   of divergence and curl operators in msu22 and psu22.                        
c                                                                               
c program history log:                                                          
c   88-11-02  joseph sela                                                       
c                                                                               
c usage:    call gozrmt (qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgtl)              
c   input argument list:                                                        
c     qlnt     - doubled scalar triangular                                      
c                array of associated legendre functions at                      
c                a given latitude.                                              
c                on input, values of qlnt are a subset of qlnv.                 
c     qlnv     - doubled vector triangular                                      
c                array of associated legendre functions at                      
c                a given latitude.                                              
c     eps      - array of function of wave number computed in epslon.           
c                eps is used only during first call to gozrmt.                  
c     lat      - latitude index.                                                
c     rcs2     - array of constants computed in glats (1/sin(lat)**2).          
c     wgtl     - weight at gaussian latitude.                                   
c                multiplier of output arrays qlnt, qdert, qlnwct.               
c                                                                               
c   output argument list:                                                       
c     qlnt     - doubled scalar triangular                                      
c                array of qlnt*n*(n+1)*1/a**2   times wgtl.                     
c     qdert    - doubled scalar triangular                                      
c                array of legendre derivatives  times wgtl.                     
c     qlnwct   - doubled scalar triangular                                      
c                array of qlnt*l*rcs2(lat)/a    times wgtl.                     
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
cra       save                                                                 
cc                                                                              
      dimension          qlnt(lnt2_)                                            
      dimension          qlnv(lnut2_)                                           
      dimension         qdert(lnt2_)                                            
      dimension           eps(jcap1_,jcap2_)                                    
      dimension        qlnwct(lnt2_)                                            
      dimension          rcs2(latg2_)                                           
cc                                                                              
      dimension           dxa(lnt2_)                                            
      dimension           dxb(lnt2_)                                            
      dimension           dxc(lnt2_)                                            
      dimension           dxd(lnt2_)                                            
cc                                                                              
      dimension         dxint(lnt2_)                                            
      dimension            dx(twoj1_,jcap2_)                                    
      dimension          deps(twoj1_,jcap2_)                                    
c     common /gozcom/ifir,jfir,dxa,dxb,dxc,dxd                                  
      common /gozcom/dxa,dxb,dxc,dxd,jfir                                       
cc                                                                              
      data ifir /0/                                                             
      if  (ifir .eq. 1)  go to 500                                              
           ifir = 1                                                             
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
             deps(2*ll-1,i) = eps(ll,i)                                         
             deps(2*ll  ,i) = eps(ll,i)                                         
  260 continue                                                                  
  280 continue                                                                  
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
      lp  = 0                                                                   
      len = twoj1_                                                              
      do 440  i=1,jcap1_                                                        
      do 420 ll=1,len                                                           
             dxc(ll+lp) = dx(ll,1)                                              
             dxd(ll+lp) = dx(ll,i) * dx(ll,i+1)                                 
  420 continue                                                                  
      lp  = lp  + len                                                           
      len = len - 2                                                             
  440 continue                                                                  
  500 continue                                                                  
cc                                                                              
cc    compute pln derivatives                                                   
      wcsa=rcs2(lat)*wgtl/rerth_                                                
      raa=wgtl/( rerth_ * rerth_ )                                              
cc                                                                              
      lp0 = 0                                                                   
      lp1 = twoj1_                                                              
      len = twoj1_ - 2                                                          
      do 640  i=1,jcap_                                                         
      do 620 ll=1,len                                                           
             qdert(ll+lp1) = qlnt(ll+lp0) * dxa(ll+lp1)                         
  620 continue                                                                  
      lp0 = lp1                                                                 
      lp1 = lp1 + len                                                           
      len = len - 2                                                             
  640 continue                                                                  
                lbegin = twoj1_ + 1                                             
      do 720 ll=lbegin,lnt2_                                                    
             qdert(ll) = qdert(ll) + qlnv(ll+twoj1_) * dxb(ll)                  
  720 continue                                                                  
      do 740 ll=1,twoj1_                                                        
             qdert(ll) = qlnv(ll+twoj1_) * dxb(ll)                              
  740 continue                                                                  
      do 760 ll=1,lnt2_                                                         
             qdert(ll) = qdert(ll) * wcsa                                       
  760 continue                                                                  
cc                                                                              
      do 820 ll=1,lnt2_                                                         
             qlnwct(ll) = qlnt(ll) * dxc(ll) * wcsa                             
  820 continue                                                                  
cc                                                                              
      do 840 ll=1,lnt2_                                                         
             qlnt(ll) = qlnt(ll) * dxd(ll) * raa                                
  840 continue                                                                  
cc                                                                              
      return                                                                    
      end                                                                       
