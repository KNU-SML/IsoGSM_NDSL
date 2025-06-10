#define JOSGB josgb
      subroutine glat(jh,slat,clat,wlat)                                        
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    glat        compute gaussian latitude functions                
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes sines of gaussian latitude by iteration.                   
c           the cosines of gaussian latitude and gaussian weights               
c           are also computed.                                                  
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call glat(jh,slat,clat,wlat)                                        
c                                                                               
c   input argument list:                                                        
c     jh       - integer number of gaussian latitudes in a hemisphere           
c                                                                               
c   output argument list:                                                       
c     slat     - real (jh) sines of (positive) gaussian latitude                
c     clat     - real (jh) cosines of gaussian latitude                         
c     wlat     - real (jh) gaussian weights for the nh                          
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension slat(jh),clat(jh),wlat(jh)                                      
      parameter(pi=3.14159265358979,c=(1.-(2./pi)**2)*0.25,eps=1.e-14)          
      parameter(jbz=50)                                                         
      dimension bz(jbz)                                                         
      parameter(jjh=(jo_+1)/2)                                                  
      dimension pk(jjh),pkm1(jjh)                                               
      data bz        / 2.4048255577,  5.5200781103,                             
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,              
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,              
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,              
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,              
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,              
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,              
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,              
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,              
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,              
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,              
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,              
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /             
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  estimate latitudes using bessel function                                     
      r=1./sqrt((2*jh+0.5)**2+c)                                                
      do j=1,min(jh,jbz)                                                        
        slat(j)=cos(bz(j)*r)                                                    
      enddo                                                                     
      do j=jbz+1,jh                                                             
        slat(j)=cos((bz(jbz)+(j-jbz)*pi)*r)                                     
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  converge until all sines of gaussian latitude are within eps                 
      spmax=1.                                                                  
      do while(spmax.gt.eps)                                                    
        spmax=0.                                                                
        do j=1,jh                                                               
          pkm1(j)=1.                                                            
          pk(j)=slat(j)                                                         
        enddo                                                                   
        do n=2,2*jh                                                             
          do j=1,jh                                                             
            pkm2=pkm1(j)                                                        
            pkm1(j)=pk(j)                                                       
            pk(j)=((2*n-1)*slat(j)*pkm1(j)-(n-1)*pkm2)/n                        
          enddo                                                                 
        enddo                                                                   
        do j=1,jh                                                               
          sp=pk(j)*(1.-slat(j)**2)/(2*jh*(pkm1(j)-slat(j)*pk(j)))               
          slat(j)=slat(j)-sp                                                    
          spmax=max(spmax,abs(sp))                                              
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute cosines and gaussian weights                                         
      do j=1,jh                                                                 
        clat(j)=sqrt(1.-slat(j)**2)                                             
        wlat(j)=2.*(1.-slat(j)**2)/(2*jh*pkm1(j))**2                            
      enddo                                                                     
      return                                                                    
      end                                                                       
