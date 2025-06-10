#define JCAP jcap
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
      subroutine sicdif(d,t,q,x,y,z,u,v,lnt2,lnoffset)
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    sicdif      semi-implicit time integration.                    
c   prgmmr: joseph sela      org: w/nmc23    date: 89-03-15                     
c                                                                               
c abstract: integrates divergence, temperature and log surface pressure         
c           semi-implicitly in time.                                            
c                                                                               
c program history log:                                                          
c   89-03-15  joseph sela                                                       
c   93-03-15  mark iredell   linear matrices passed in common                   
c   99-03-15  hann-ming henry juang  pass offset and length
c                                                                               
c usage:    call sicdif(d,t,q,x,y,z,u,v)                                        
c   input argument list:                                                        
c     d        - divergence at time t-dt                                        
c     t        - temperature at time t-dt                                       
c     q        - ln(psfc) at time t-dt                                          
c     x        - divergence nonlinear tendency at time t                        
c     y        - temperature nonlinear tendency at time t                       
c     z        - ln(psfc) nonlinear tendency at time t                          
c                                                                               
c   output argument list:                                                       
c     x        - divergence at time t+dt                                        
c     y        - temperature at time t+dt                                       
c     z        - ln(psfc) at time t+dt                                          
c     u        - work array                                                     
c     v        - work array                                                     
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
#include <paramodel.h>
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
      parameter(km=levs_,jcap=jcap_,lnt22=LNT22S)                    
      real d(lnt22,km),t(lnt22,km),q(lnt22)                                     
      real x(lnt22,km),y(lnt22,km),z(lnt22),u(lnt22,km),v(lnt22,km)             
#include <comio.h>
      common/comsic/ dt,gvdt(km),svdt(km),amdt(km,km),bmdt(km,km),              
     1               dm(km,km,0:jcap)                                           
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  explicitly integrate lnps and temperature halfway in time.                   
      do i=1,lnt2                                                               
        z(i)=q(i)+dt*z(i)                                                       
      enddo                                                                     
      do k=1,km                                                                 
        do i=1,lnt2                                                             
          y(i,k)=t(i,k)+dt*y(i,k)                                               
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute linear dependence of divergence on lnps and temperature.             
c  explicitly integrate divergence halfway including linear terms.              
#ifdef CRAY_THREAD
cmic$ do all autoscope                                                          
#endif
      do k=1,km                                                                 
        do i=1,lnt2                                                             
          v(i,k)=0.                                                             
        enddo                                                                   
#ifdef CRAY_THREAD
cfpp$ unroll l                                                                  
#endif
        do j=1,km                                                               
          do i=1,lnt2                                                           
            v(i,k)=v(i,k)+amdt(k,j)*y(i,j)                                      
          enddo                                                                 
        enddo                                                                   
        do i=1,lnt2                                                             
          ii=lnoffset+i
          u(i,k)=d(i,k)+dt*x(i,k)+snnp1(ii)*(v(i,k)+gvdt(k)*z(i))
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  solve helmholz equation for semi-implicit divergence.                        
#ifdef CRAY_THREAD
cmic$ do all autoscope                                                          
#endif
      do k=1,km                                                                 
        do i=1,lnt2                                                             
          v(i,k)=0.                                                             
        enddo                                                                   
#ifdef CRAY_THREAD
cfpp$ unroll l                                                                  
#endif
        do j=1,km                                                               
          do i=1,lnt2,2                                                         
            ii=lnoffset+i
            n=ndex(ii)
            v(i,k)=v(i,k)+dm(k,j,n)*u(i,j)                                      
            v(i+1,k)=v(i+1,k)+dm(k,j,n)*u(i+1,j)                                
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  back solve for lnps.                                                         
#ifdef CRAY_THREAD
cfpp$ unroll l                                                                  
#endif
      do j=1,km                                                                 
        do i=1,lnt2                                                             
          z(i)=z(i)+svdt(j)*v(i,j)                                              
        enddo                                                                   
      enddo                                                                     
      do i=1,lnt2                                                               
        z(i)=2*z(i)-q(i)                                                        
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  back solve for temperature and divergence.                                   
#ifdef CRAY_THREAD
cmic$ do all autoscope                                                          
#endif
      do k=1,km                                                                 
#ifdef CRAY_THREAD
cfpp$ unroll l                                                                  
#endif
        do j=1,km                                                               
          do i=1,lnt2                                                           
            y(i,k)=y(i,k)+bmdt(k,j)*v(i,j)                                      
          enddo                                                                 
        enddo                                                                   
        do i=1,lnt2                                                             
          y(i,k)=2*y(i,k)-t(i,k)                                                
          x(i,k)=2*v(i,k)-d(i,k)                                                
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
