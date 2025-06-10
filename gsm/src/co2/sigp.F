#define LEVS levs
      subroutine sigp(pstar,pd,gtemp,t41,t42,t43,t44,sglvnu,siglnu,             
     1                lread)                                                    
#include <paramodel.h>
      parameter (kd=levs_)                                                      
      parameter (kp=kd+1,km=kd-1,kp2=kd+2)                                      
      dimension q(kd),qmh(kp),pd(kp2),plm(kp),gtemp(kp),pdt(kp2)                
      dimension ci(kp),sglvnu(kp),del(kd),siglnu(kd),cl(kd),rpi(km)             
      dimension idate(4)                                                        
      dimension t41(kp2,2),t42(kp),                                             
     1          t43(kp2,2),t44(kp)                                              
ccc   18 level sigmas for nmc mrf(new) model                                    
ccc   data q/.021,.074,.124,.175,.225,.275,.325,.375,.425,.497,                 
ccc  1       .594,.688,.777,.856,.920,.960,.981,.995/                           
c     for sigma models,q=sigma,qmh=0.5(q(i)+q(i+1),                             
c     pd=q*pss,plm=qmh*pss.pss=surface pressure(spec.)                          
c                                                                               
c.....   get nmc sigma structure                                                
      call newsig(ci,sglvnu,del,siglnu,cl,rpi)                                  
      do 77 k=1,kd                                                              
         q(k) = siglnu(kd+1-k)                                                  
   77 continue                                                                  
      pss=    1013250.                                                          
      qmh(1)=0.                                                                 
      qmh(kp)=1.                                                                
      do 1 k=2,kd                                                               
      qmh(k)=0.5*(q(k-1)+q(k))                                                  
1     continue                                                                  
      pd(1)=0.                                                                  
      pd(kp2)=pss                                                               
      do 2 k=2,kp                                                               
      pd(k)=q(k-1)*pss                                                          
2     continue                                                                  
      plm(1)=0.                                                                 
      do 3 k=1,km                                                               
      plm(k+1)=0.5*(pd(k+1)+pd(k+2))                                            
3     continue                                                                  
      plm(kp)=pss                                                               
      do 4 k=1,kd                                                               
      gtemp(k)=pd(k+1)**0.2*(1.+pd(k+1)/30000.)**0.8/1013250.                   
4     continue                                                                  
      gtemp(kp)=0.                                                              
      write (6,100) (gtemp(k),k=1,kd)                                           
      write (6,100) (pd(k),k=1,kp2)                                             
      write (6,100) (plm(k),k=1,kp)                                             
c***tapes 41,42 are output to the co2 interpolation program (ps=1013mb)         
c  the following puts p-data into mb                                            
      do 11 i=1,kp                                                              
      pd(i)=pd(i)*1.0e-3                                                        
      plm(i)=plm(i)*1.0e-3                                                      
11    continue                                                                  
      pd(kp2)=pd(kp2)*1.0e-3                                                    
ccc         write (41,101) (pd(k),k=1,kp2)                                      
ccc         write (41,101) (plm(k),k=1,kp)                                      
ccc         write (42,101) (plm(k),k=1,kp)                                      
      do 300 k=1,kp2                                                            
       t41(k,1) = pd(k)                                                         
  300 continue                                                                  
      do 301 k=1,kp                                                             
       t41(k,2) = plm(k)                                                        
       t42(k) = plm(k)                                                          
  301 continue                                                                  
c***store as pdt,so that right pd is returned to ptz                            
      do 12 i=1,kp2                                                             
      pdt(i)=pd(i)                                                              
12    continue                                                                  
c***second pass: pss=810mb,gtemp not computed                                   
      pss=0.8*1013250.                                                          
      qmh(1)=0.                                                                 
      qmh(kp)=1.                                                                
      do 201 k=2,kd                                                             
      qmh(k)=0.5*(q(k-1)+q(k))                                                  
201   continue                                                                  
      pd(1)=0.                                                                  
      pd(kp2)=pss                                                               
      do 202 k=2,kp                                                             
      pd(k)=q(k-1)*pss                                                          
202   continue                                                                  
      plm(1)=0.                                                                 
      do 203 k=1,km                                                             
      plm(k+1)=0.5*(pd(k+1)+pd(k+2))                                            
203   continue                                                                  
      plm(kp)=pss                                                               
      write (6,100) (pd(k),k=1,kp2)                                             
      write (6,100) (plm(k),k=1,kp)                                             
c***tapes 43,44 are output to the co2 interpolation program(ps=810 mb)          
c  the following puts p-data into mb                                            
      do 211 i=1,kp                                                             
      pd(i)=pd(i)*1.0e-3                                                        
      plm(i)=plm(i)*1.0e-3                                                      
211   continue                                                                  
      pd(kp2)=pd(kp2)*1.0e-3                                                    
ccc       write (43,101) (pd(k),k=1,kp2)                                        
ccc       write (43,101) (plm(k),k=1,kp)                                        
ccc       write (44,101) (plm(k),k=1,kp)                                        
      do 302 k=1,kp2                                                            
       t43(k,1) = pd(k)                                                         
  302 continue                                                                  
      do 303 k=1,kp                                                             
       t43(k,2) = plm(k)                                                        
       t44(k) = plm(k)                                                          
  303 continue                                                                  
c***restore pd                                                                  
      do 212 i=1,kp2                                                            
      pd(i)=pdt(i)                                                              
212   continue                                                                  
100   format (1x,5e20.13)                                                       
101   format (5e16.9)                                                           
      return                                                                    
      end                                                                       
