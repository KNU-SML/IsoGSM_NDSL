#define JCAP jcap
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
      subroutine impadj_hyb(d,t,q,x,y,z,u,v,lnt2,lnoffset)
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    impadj      implicit adjustment of physics tendencies.         
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: extends the semi-implicit time integration to include               
c           the phyical forcing terms computed in gloopb.                       
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell                                                      
c   93-03-15  mark iredell   change argument list                               
c   93-06-22  hann-ming henry juang  pass offset and length for mpp
c                                                                               
c usage:    call impadj(d,t,q,x,y,z,u,v)                                        
c   input argument list:                                                        
c     d        - divergence before adjustment                                   
c     t        - temperature before adjustment                                  
c     q        - ln(psfc) before adjustment                                     
c     x        - divergence tendency adjustment                                 
c     y        - temperature tendency adjustment                                
c     z        - ln(psfc) tendency adjustment                                   
c                                                                               
c   output argument list:                                                       
c     d        - divergence adjusted                                            
c     t        - temperature adjusted                                           
c     q        - ln(psfc) adjusted                                              
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
c  compute linear dependence of divergence on lnps and temperature.             
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
          u(i,k)=x(i,k)+snnp1(ii)*(v(i,k)+gvdt(k)*z(i))
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
c          q(i)=q(i)+svdt(j)*v(i,j)                                              
          q(i)=q(i)-svdt(j)*v(i,j)                                              
        enddo                                                                   
      enddo                                                                     
      do i=1,lnt2                                                               
        q(i)=q(i)+z(i)                                                          
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
c            t(i,k)=t(i,k)+bmdt(k,j)*v(i,j)                                      
            t(i,k)=t(i,k)-bmdt(k,j)*v(i,j)                                      
          enddo                                                                 
        enddo                                                                   
        do i=1,lnt2                                                             
          t(i,k)=t(i,k)+y(i,k)                                                  
          d(i,k)=d(i,k)+v(i,k)                                                  
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
