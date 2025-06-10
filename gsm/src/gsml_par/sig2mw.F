#define IGRD igrd
#define IGRD1 igrd1
#define LEVS levs
#define NPES npes
      subroutine sig2mw(im,ix,km,sl,ps,u,v,t,pmw,umw,vmw,tmw)                   
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    sig2mw      sigma to maxwind interpolation                    
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: locates the maximum wind speed level (maxwind level) and            
c   returns the wind speed, components and pressure at that level.              
c   the maxwind level is restricted to be between 50kpa and 7kpa.               
c   the maxwind level is identified by cubic spline interpolation               
c   of the wind speeds in log pressure.                                         
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c   94-08-01  iredell   return maxwind temperature too                          
c                                                                               
c usage:    call sig2mw(im,ix,km,sl,ps,u,v,t,pmw,umw,vmw,tmw)                   
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
c     pmw      - real (im) maxwind pressure in kpa                              
c     umw      - real (im) maxwind zonal wind in m/s                            
c     vmw      - real (im) maxwind merid wind in m/s                            
c     tmw      - real (im) maxwind temperature in k                             
c                                                                               
c subprograms called:                                                           
c   spcoef       compute second derivatives for cubic spline                    
c   spfmax       determine maximum value of cubic spline                        
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension sl(im,km),ps(im)                                                   
      dimension u(ix,km),v(ix,km),t(ix,km)                                      
      dimension pmw(im),umw(im),vmw(im),tmw(im)                                 
      dimension spdmw(im),smw(im),s(im,levs_)                              
      dimension spd(im,levs_),d2spd(im,levs_)                           
      parameter(pmwbot=500.e-1,pmwtop=70.e-1)                                   
c
c  fix vertical coordinate proportional to log pressure                         
c  and calculate wind speeds between pmwbot and pmwtop                          
c
      do k=1,km                                                                 
        do i=1,im                                                               
          s(i,k)=-log(sl(i,k))                                                        
          p=sl(i,k)*ps(i)                                                         
          if(p.le.pmwbot.and.p.ge.pmwtop) then                                  
            spd(i,k)=sqrt(u(i,k)**2+v(i,k)**2)                                  
          else                                                                  
            spd(i,k)=0.                                                         
          endif                                                                 
        enddo                                                                   
      enddo                                                                     
c
c  use spline routines to determine maxwind level and wind speed                
c
      call spcoef(im,km,s,spd,d2spd)                                            
      call spfmax(im,km,s,spd,d2spd,smw,pmw,spdmw)                              
c
c  compute maxwind pressure and wind components                                 
c
      do i=1,im                                                                 
        pmw(i)=exp(-pmw(i))*ps(i)                                               
        k=int(smw(i))                                                           
        if(float(k).eq.smw(i)) then                                             
          ub=u(i,k)                                                             
          vb=v(i,k)                                                             
          tb=t(i,k)                                                             
        else                                                                    
          ub=(k+1-smw(i))*u(i,k)+(smw(i)-k)*u(i,k+1)                            
          vb=(k+1-smw(i))*v(i,k)+(smw(i)-k)*v(i,k+1)                            
          tb=(k+1-smw(i))*t(i,k)+(smw(i)-k)*t(i,k+1)                            
        endif                                                                   
        spdb=sqrt(ub**2+vb**2)                                                  
        umw(i)=ub*spdmw(i)/spdb                                                 
        vmw(i)=vb*spdmw(i)/spdb                                                 
        tmw(i)=tb                                                               
      enddo                                                                     
c
      return                                                                    
      end                                                                       
