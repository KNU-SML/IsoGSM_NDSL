#define LEVS levs
      subroutine sig2tp(im,ix,km,sl,                                            
     &                  ps,u,v,t,                                               
     &                  ptp,utp,vtp,ttp,shtp,stp)                               
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    sig2tp      sigma to tropopause interpolation                 
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: locates the tropopause pressure level and interpolates              
c   the winds and temperature and wind shear to the tropopause.                 
c   the tropopause is identified by the lowest level above 450 mb               
c   where the temperature lapse rate -dt/dz becomes less than 2 k/km.           
c   the tropopause is not allowed higher than 85 mb.                            
c   interpolations are done linearly in log of pressure.                        
c                                                                               
c program history log:                                                          
c   92-10-31  mccalla,iredell                                                   
c                                                                               
c usage:    call sig2tp(im,ix,km,sl,                                            
c    &                  ps,u,v,t,                                               
c    &                  ptp,utp,vtp,ttp,shtp,stp)                               
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of sigma levels                                 
c     sl       - real (km) sigma values                                         
c     ps       - real (im) surface pressure in kpa                              
c     u        - real (ix,km) zonal wind in m/s                                 
c     v        - real (ix,km) merid wind in m/s                                 
c     t        - real (ix,km) temperature in k                                  
c                                                                               
c   output argument list:                                                       
c     ptp      - real (im) tropopause pressure in kpa                           
c     utp      - real (im) tropopause zonal wind in m/s                         
c     vtp      - real (im) tropopause merid wind in m/s                         
c     ttp      - real (im) tropopause temperature in k                          
c     shtp     - real (im) tropopause wind speed shear in (m/s)/m               
c     stp      - real (im) tropopause sigma layer number                        
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension sl(km),ps(im)                                                   
      dimension u(ix,km),v(ix,km),t(ix,km)                                      
      dimension ptp(im),utp(im),vtp(im),ttp(im),shtp(im),stp(im)                
      dimension asl(levs_)                                                      
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 )                                   
      parameter(rog=rd/g)                                                       
      parameter(ptbot=450.e-1,pttop=85.e-1,gamt=2.e-3)                          
      fgamma(k)=(t(i,k-1)-t(i,k+1))/(rog*t(i,k)*(asl(k-1)-asl(k+1)))            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  identify tropopause as first layer above ptbot but below pttop               
c  where the temperature lapse rate drops below gamt                            
cim  stp is real interpolated sigma layer number of tropopause                  
      do k=1,km                                                                 
      asl(k)=log(sl(k))                                                         
      enddo                                                                     
      do i=1,im                                                                 
        k=3                                                                     
        pu=ps(i)*sl(k)                                                          
        dowhile(k.lt.km-1.and.pu.gt.ptbot)                                      
          k=k+1                                                                 
          pu=ps(i)*sl(k)                                                        
        enddo                                                                   
        gamd=fgamma(k-1)                                                        
        gamd=max(gamd,gamt)                                                     
        gamu=fgamma(k)                                                          
        dowhile(k.lt.km-1.and.pu.gt.pttop.and.gamu.gt.gamt)                     
          k=k+1                                                                 
          pu=ps(i)*sl(k)                                                        
          gamd=gamu                                                             
          gamu=fgamma(k)                                                        
        enddo                                                                   
        gamu=min(gamu,gamt)                                                     
        stp(i)=k-(gamt-gamu)/(gamd-gamu)                                        
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  interpolate tropopause pressure, temperature, winds and wind shear           
c  tropopause pressure is constrained to be between ptbot and pttop             
      do i=1,im                                                                 
        kd=stp(i)                                                               
        ku=kd+1                                                                 
        wu=stp(i)-kd                                                            
        dlp=asl(ku)-asl(kd)                                                     
        ptp(i)=ps(i)*sl(kd)*exp(wu*dlp)                                         
        if(ptp(i).gt.ptbot) then                                                
          wu=wu+log(ptbot/ptp(i))/dlp                                           
          ptp(i)=ptbot                                                          
        elseif(ptp(i).lt.pttop) then                                            
          wu=wu+log(pttop/ptp(i))/dlp                                           
          ptp(i)=pttop                                                          
        endif                                                                   
        ttp(i)=t(i,kd)+wu*(t(i,ku)-t(i,kd))                                     
        utp(i)=u(i,kd)+wu*(u(i,ku)-u(i,kd))                                     
        vtp(i)=v(i,kd)+wu*(v(i,ku)-v(i,kd))                                     
        spdd=sqrt(u(i,kd)**2+v(i,kd)**2)                                        
        spdu=sqrt(u(i,ku)**2+v(i,ku)**2)                                        
        shtp(i)=(spdu-spdd)/(rog*0.5*(t(i,ku)+t(i,kd))*dlp)                     
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
