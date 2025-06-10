#define IO io
#define KO ko
#define LEVS levs
#include <define.h>
      subroutine sig2p(im,ix,km,si,sl,                                          
     &                 ps,us,vs,os,zs,zi,ts,rs,qs,cls,prs,o3s,                              
     &                 ko,po,lpcl,lppr,lpo3,up,vp,op,zp,tp,rp,
     &                 qp,clp,prp,o3p)                                 
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    sig2p       sigma to pressure interpolation                   
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: interpolates winds, omega, height, temperature and humidity         
c   from the sigma coordinate system to the mandatory pressure levels.          
c   assumes that relative humidity, temperature, geopotential heights,          
c   wind components and vertical velocity vary linearly in the vertical         
c   with the log of pressure.  underground heights are obtained using           
c   the shuell method and underground temperatures are obtained using           
c   a constant moist adiabatic lapse rate.  heights above the top sigma         
c   level are integrated hydrostatically.  otherwise fields are held            
c   constant outside the sigma structure and no extrapolation is done.          
c                                                                               
c program history log:                                                          
c   92-10-31  sela,newell,gerrity,ballish,deaven,iredell                        
c                                                                               
c usage:    call sig2p(im,ix,km,si,sl,                                          
c    &                 ps,us,vs,os,zs,zi,ts,rs,qs,                              
c    &                 ko,po,up,vp,op,zp,tp,rp,sp)                              
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of sigma levels                                 
c     si       - real (km+1) sigma interface values                             
c     sl       - real (km) sigma values                                         
c     ps       - real (im) surface pressure in kpa                              
c     us       - real (ix,km) zonal wind in m/s                                 
c     vs       - real (ix,km) merid wind in m/s                                 
c     os       - real (ix,km) vertical velocity in pa/s                         
c     zs       - real (ix,km) heights on the full levels in m                   
c     zi       - real (ix,km) heights on the interfaces in m                    
c     ts       - real (ix,km) temperature in k                                  
c     rs       - real (ix,km) relative humidity in percent                      
c     qs       - real (ix,km) specific humidity in kg/kg                        
c     ko       - integer number of pressure levels                              
c     po       - real (ko) mandatory pressures in kpa                           
c                                                                               
c   output argument list:                                                       
c     up       - real (ix,ko) zonal wind in m/s                                 
c     vp       - real (ix,ko) merid wind in m/s                                 
c     op       - real (ix,ko) vertical velocity in pa/s                         
c     zp       - real (ix,ko) heights in m                                      
c     tp       - real (ix,ko) temperature in k                                  
c     rp       - real (ix,ko) relative humidity in percent                      
c                                                                               
c subprograms called:                                                           
c   isrchfltx - find first value in an array less than target value              
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension si(ix,km),sl(ix,km),ps(im)                                            
      dimension us(ix,km),vs(ix,km),os(ix,km)                                   
      dimension zs(ix,km),zi(ix,km),ts(ix,km),rs(ix,km),qs(ix,km)               
      dimension cls(ix,km),prs(ix,km),o3s(ix,km)
      dimension po(ko)                                                          
      logical lpcl,lppr,lpo3
      dimension up(ix,ko),vp(ix,ko),op(ix,ko)                                   
      dimension zp(ix,ko),tp(ix,ko),rp(ix,ko),qp(ix,ko)
      dimension clp(ix,ko),prp(ix,ko),o3p(ix,ko)
      dimension sp(2*io_+6,ko_)                                                 
      dimension asi(levs_),asl(levs_),apo(ko_),aps(2*io_)                       
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 ,rv= 4.6150e+2 )                    
      parameter(rog=rd/g,fvirt=rv/rd-1.)                                        
      parameter(gammam=-6.5e-3,zshul=75.,tvshul=290.66)                         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute log pressures for interpolation                                      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  determine sigma layers bracketing pressure layer                             
c  and interpolate to obtain real sigma layer number                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  interpolate sigma to pressure                                                
      do i=1,im
        do k=1,km
          asi(k)=log(si(i,k))
          asl(k)=log(sl(i,k))
        enddo
        do k=1,ko
          apo(k)=log(po(k))
        enddo
        aps(i)=log(ps(i))

        kd=1
        do k=1,ko
          ask=apo(k)-aps(i)
          kd=kd+isrchfltx(km-kd-1,asl(kd+1),1,ask)-1
          sp(i,k)=kd+(asl(kd)-ask)/(asl(kd)-asl(kd+1))
        enddo
        do k=1,ko
          ask=apo(k)-aps(i)
c  below ground use shuell method to obtain height, constant lapse rate         
c  to obtain temperature, and hold other fields constant                        
          if(ask.gt.0.) then                                                    
            up(i,k)=us(i,1)                                                     
            vp(i,k)=vs(i,1)                                                     
            op(i,k)=os(i,1)                                                     
            tvsf=ts(i,1)*(1.+fvirt*qs(i,1))-gammam*(zs(i,1)-zi(i,1))            
            if(zi(i,1).gt.zshul) then                                           
              tvsl=tvsf-gammam*zi(i,1)                                          
              if(tvsl.gt.tvshul) then                                           
                if(tvsf.gt.tvshul) then                                         
                  tvsl=tvshul-5.e-3*(tvsf-tvshul)**2                            
                else                                                            
                  tvsl=tvshul                                                   
                endif                                                           
              endif                                                             
              gammas=(tvsf-tvsl)/zi(i,1)                                        
            else                                                                
              gammas=0.                                                         
            endif                                                               
            part=rog*ask                                                        
            zp(i,k)=zi(i,1)-tvsf*part/(1.+0.5*gammas*part)                      
            tp(i,k)=ts(i,1)+gammam*(zp(i,k)-zs(i,1))                            
            rp(i,k)=rs(i,1)                                                     
            qp(i,k)=qs(i,1)                                                     
            if(lpcl) clp(i,k)=cls(i,1)
            if(lppr) prp(i,k)=prs(i,1)
            if(lpo3) o3p(i,k)=o3s(i,1)
c  above top sigma ground integrate height hydrostatically                      
c  and hold other fields constant                                               
          elseif(sp(i,k).ge.km) then                                            
            up(i,k)=us(i,km)                                                    
            vp(i,k)=vs(i,km)                                                    
            op(i,k)=os(i,km)                                                    
            tvkm=ts(i,km)*(1.+fvirt*qs(i,km))                                   
            zp(i,k)=zs(i,km)+rog*tvkm*(asl(km)-ask)                             
            tp(i,k)=ts(i,km)                                                    
            rp(i,k)=rs(i,km)                                                    
            qp(i,k)=qs(i,km)                                                    
            if(lpcl) clp(i,k)=cls(i,km)
            if(lppr) prp(i,k)=prs(i,km)
            if(lpo3) o3p(i,k)=o3s(i,km)
c  within sigma structure, interpolate fields linearly in log pressure          
c  between bracketing full sigma layers except heights are interpolated         
c  between the nearest full sigma layer and the nearest sigma interface         
          else                                                                  
            kd=max(int(sp(i,k)),1)                                              
            ku=kd+1                                                             
            wu=sp(i,k)-kd                                                       
            wd=1.-wu                                                            
            up(i,k)=wu*us(i,ku)+wd*us(i,kd)                                     
            vp(i,k)=wu*vs(i,ku)+wd*vs(i,kd)                                     
            op(i,k)=wu*os(i,ku)+wd*os(i,kd)                                     
            ki=int(sp(i,k))+1                                                   
            di=asi(ki)-ask                                                      
            kl=nint(ki-0.5+sign(0.5,di))                                        
            wl=di/(asi(ki)-asl(kl))                                             
            wi=1.-wl                                                            
            zp(i,k)=wi*zi(i,ki)+wl*zs(i,kl)                                     
            tp(i,k)=wu*ts(i,ku)+wd*ts(i,kd)                                     
            rp(i,k)=wu*rs(i,ku)+wd*rs(i,kd)                                     
            qp(i,k)=wu*qs(i,ku)+wd*qs(i,kd)
            if(lpcl) clp(i,k)=wu*cls(i,ku)+wd*cls(i,kd)
            if(lppr) prp(i,k)=wu*prs(i,ku)+wd*prs(i,kd)
            if(lpo3) o3p(i,k)=wu*o3s(i,ku)+wd*o3s(i,kd)
          endif                                                                 
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
