#undef sun
#define IO io
#include <define.h>
      subroutine sundry(im,ix,km,nsun,kslp,si,zs,ps,q,zm,sun)                   
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    sundry      compute sundry fields                              
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: computes sundry fields from sigma level winds, omega,               
c   temperature, moisture, and relative humidity.                               
c   relative humidity is perversely averaged across sigma layers.               
c   the current nsun=4 sundry fields are:                                      
c     1) surface pressure                                                       
c     2) precipitable water                                                     
c     3) surface orography                                                      
c     4) sea level pressure                                                     
c                                                                               
c subprograms called:                                                           
c   sig2tp       interpolate sigma to tropopause level                          
c   sig2mw       interpolate sigma to maxwind level                             
c   liftix       compute best lifted index                                      
c                                                                               
c program history log:                                                          
c   92-10-31  mccalla,iredell                                                   
c   96-01-18  kanamistu  greatly simplified                                     
c                                                                               
c usage:    call sundry(im,ix,km,kslp,rnhour,si,                                
c    &                  zs,ps,q,zm,sun)                                         
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of levels                                       
c     kslp     - integer (2) relevant pressure levels for slp                   
c     rnhour   - real hours rain accumulated                                    
c     si       - real (km) sigma interfaces                                     
c     zs       - real (im) surface orography in m                               
c     ps       - real (im) surface pressure in kpa                              
c     q        - real (ix,km) specific humidity in kg/kg                        
c     zm       - real (ix,*) height on pressure surface in m                    
c                                                                               
c   output argument list:                                                       
c     sun      - real (ix,nsun) sundry fields given above                       
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension kslp(2)                                                         
      dimension si(ix,km+1)                                                        
      dimension zs(im),ps(im)                                                   
c                                                                               
      dimension q(ix,km)                                                        
      dimension zm(ix,km)                                                       
      dimension sun(ix,nsun)                                                    
      dimension sk(2*io_),dum(2*io_)                                            
      parameter(g= 9.8000e+0 ,cp= 1.0046e+3 )                                   
      parameter(rd= 2.8705e+2 ,rv= 4.6150e+2 )                                  
      parameter(pm1=1.e5,tm1=287.45,zm1=113.,zm2=5572.)                         
      parameter(rk=rd/cp,fslp=g*(zm2-zm1)/(rd*tm1))                             
c                                                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  surface pressure                                                             
#ifdef DBG
      call maxmin(ps,im,1,1,1,'ps')                                               
#endif
      do i=1,im                                                                 
      sun(i,1)=ps(i)*1.e3                                                       
      enddo                                                                     
#ifdef DBG
      call maxmin(sun(1,1),ix,1,1,1,'sun1')                                         
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  precipitable water                                                           
#ifdef DBG
      print *,'q'                                                               
      call maxmin(q(1,1),ix,km,1,1,'q')                                         
#endif
      do i=1,im                                                                 
      sun(i,2)=0.                                                               
      enddo                                                                     
      do k=1,km                                                                 
        do i=1,im                                                               
          ds=si(i,k)-si(i,k+1)                                                        
          sun(i,2)=sun(i,2)+q(i,k)*ds                                           
        enddo                                                                   
      enddo                                                                     
      do i=1,im                                                                 
      sun(i,2)=sun(i,2)*ps(i)*1.e3/g                                            
      enddo                                                                     
#ifdef DBG
      print *,'sun2'                                                            
      call maxmin(sun(1,2),ix,1,1,1,'sun2')                                         
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  surface orography                                                            
#ifdef DBG
      print *,'zs'                                                              
      call maxmin(zs,im,1,1,1,'zs')                                               
#endif
      do i=1,im                                                                 
      sun(i,3)=zs(i)                                                            
      enddo                                                                     
#ifdef DBG
      print *,'sun3'                                                            
      call maxmin(sun(1,3),ix,1,1,1,'sun3')                                         
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  sea level pressure                                                           
      if(kslp(1).gt.0.and.kslp(2).gt.0) then                                    
      k1=kslp(1)                                                                
      k2=kslp(2)                                                                
#ifdef DBG
      print *,'k1,k2=',k1,k2                                                    
      print *,'zm(k1) k1=',k1                                                   
      call maxmin(zm(1,k1),ix,1,1,1,'zm1')
      print *,'zm(k2) k2=',k2                                                   
      call maxmin(zm(1,k2),ix,1,1,1,'zm2')                                         
#endif
        do i=1,im                                                               
        sun(i,4)=pm1*exp(fslp*zm(i,k1)/(zm(i,k2)-zm(i,k1)))                     
        enddo                                                                   
      else                                                                      
        do i=1,im                                                               
        sun(i,4)=0.                                                             
        enddo                                                                   
      endif                                                                     
#ifdef DBG
      print *,'sun4'                                                            
      call maxmin(sun(1,4),ix,1,1,1,'sun4')                                         
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
