#define IO io
#include <define.h>
      subroutine omega(im,ix,km,si,sl,                                          
#ifdef SCM
     &                 w,o)
#else
     &                 ps,psx,psy,d,u,v,o)                                      
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    omega       calculate pressure vertical velocity               
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: calculates pressure vertical velocity omega as a function           
c   of surface pressure, surface pressure gradients, and divergence             
c   and wind components on the sigma surfaces.  the formula for omega           
c   is derived from the continuity equation                                     
c     o=(sig*v.grad(lnps)-sum((d+v.grad(lnps))*dsig))*ps*1.e3                   
c   where the sum is taken from the top of the atmosphere.                      
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call omega(im,ix,km,si,sl,                                          
c    &                 ps,psx,psy,d,u,v,o)                                      
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of levels                                       
c     si       - real (km+1) sigma interface values                             
c     sl       - real (km) sigma values                                         
c     ps       - real (im) surface pressure in kpa                              
c     psx      - real (im) zonal gradient of log pressure in 1/m                
c     psy      - real (im) merid gradient of log pressure in 1/m                
c     d        - real (ix,km) divergence in 1/s                                 
c     u        - real (ix,km) zonal wind in m/s                                 
c     v        - real (ix,km) merid wind in m/s                                 
c                                                                               
c   output argument list:                                                       
c     o        - real (ix,km) pressure vertical velocity in pa/s                
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension si(ix,km+1),sl(ix,km)                                                 
#ifdef SCM
      dimension w(ix,km),o(ix,km)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do k=1,km
       do i=1,im

        o(i,k)=w(i,k)*1.e3

       enddo

      enddo
#else
      dimension ps(im),psx(im),psy(im)                                          
      dimension d(ix,km),u(ix,km),v(ix,km),o(ix,km)                             
      dimension sum(io_*2)                                                      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do i=1,im                                                                 
        sum(i)=0.                                                               
      enddo                                                                     
      do k=km,1,-1                                                              
        do i=1,im                                                               
          vgradp=u(i,k)*psx(i)+v(i,k)*psy(i)                                    
          gradpv=vgradp+d(i,k)                                                  
          sum(i)=sum(i)+gradpv*(sl(i,k)-si(i,k+1))                                  
          o(i,k)=(vgradp*sl(i,k)-sum(i))*ps(i)*1.e3                               
          sum(i)=sum(i)+gradpv*(si(i,k)-sl(i,k))                                    
        enddo                                                                   
      enddo                                                                     
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
