#define JO jo
#include <define.h>
      subroutine elat(jh,slat,clat,wlat)                                        
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    elat        compute equally-spaced latitude functions          
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes sines and cosines and gaussian weights                     
c           of equally-spaced latitudes from pole to equator.                   
c           the weights are computed based on ellsaesser (jam,1966).            
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c   93-12-28  iredell  modified weights based on ellsaesser                     
c   95-??-??  kanamitsu modified to run on workstations                         
c                                                                               
c usage:    call elat(jh,slat,clat,wlat)                                        
c                                                                               
c   input argument list:                                                        
c     jh       - integer number of latitudes in a hemisphere                    
c                                                                               
c   output argument list:                                                       
c     slat     - real (jh) sines of latitude                                    
c     clat     - real (jh) cosines of latitude                                  
c     wlat     - real (jh) gaussian weights                                     
c                                                                               
c subprograms called:                                                           
c   minv         solves full matrix problem                                     
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension slat(jh),clat(jh),wlat(jh)                                      
      parameter(jjh=(jo_+1)/2)                                                  
      dimension awork(jjh,jjh+1)                                                
#define DEFAULT
#ifdef MINV
#undef DEFAULT
      dimension bwork(jjh*2)                                                    
#endif
#ifdef DEFAULT
      dimension iwork(jjh*2)                                                    
#endif
      parameter(pi=3.14159265358979)                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      dlat=0.5*pi/(jh-1)                                                        
      slat(1)=1.                                                                
      clat(1)=0.                                                                
      do j=2,jh-1                                                               
        slat(j)=cos((j-1)*dlat)                                                 
        clat(j)=sin((j-1)*dlat)                                                 
      enddo                                                                     
      slat(jh)=0.                                                               
      clat(jh)=1.                                                               
      do js=1,jh                                                                
        do j=1,jh                                                               
          awork(js,j)=cos(2*(js-1)*(j-1)*dlat)                                  
        enddo                                                                   
      enddo                                                                     
#define DEFAULT
#ifdef MINV
#undef DEFAULT
        do js=1,jh                                                              
          awork(js,jh+1)=-1./(4*(js-1)**2-1)                                    
        enddo                                                                   
        call minv (awork,jh,jh,bwork,da,1.e-12,1,0)                             
        do j=1,jh                                                               
          wlat(j)=awork(j,jh+1)                                                 
        enddo                                                                   
#endif
#ifdef DEFAULT
      call iminv(awork,jh,det,iwork(1),iwork(jh+1))                             
      do j=1,jh                                                                 
        wlat(j)=0.                                                              
      enddo                                                                     
      do j=1,jh                                                                 
        do jj=1,jh                                                              
          wlat(j)=wlat(j)+awork(jj,j)*(-1./(4.*(jj-1)**2-1))                    
        enddo                                                                   
      enddo                                                                     
#endif
      return                                                                    
      end                                                                       
