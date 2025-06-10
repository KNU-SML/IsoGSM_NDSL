#define JCAP jcap
#define LEVS levs
#define LNT2 lnt2
#define RD rd
#define RERTH rerth

      subroutine gsicdf(deltim,am,bm,gv,sv,cm)                                  
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    gsicdf      setup for semi-implicit time integration.          
c   prgmmr: joseph sela      org: w/nmc23    date: 89-03-15                     
c                                                                               
c abstract: computes matrix inverse of rhs divergence autodependence            
c           in the semi-implicit treatment of the gravity wave modes.           
c                                                                               
c program history log:                                                          
c   93-03-15  mark iredell                                                      
c                                                                               
c usage:    call gsicdf(deltim,am,bm,gv,sv,cm)                                  
c   input argument list:                                                        
c     deltim   - timestep                                                       
c     am       - div dependence on temp (hydrostatic)                           
c     bm       - temp dependence on div (energy conversion)                     
c     gv       - div dependence on lnps (pressure gradient)                     
c     sv       - lnps dependence on div (continuity)                            
c     cm       - div autodependence (cm=sv*gv+am*bm)                            
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(km=levs_,jcap=jcap_,lnt2=lnt2_)                                 
      dimension am(km*km),bm(km*km),sv(km),gv(km),cm(km*km)                     
      parameter(rd=rd_,rerth=rerth_,raa=rd/(rerth**2))                          
      parameter(tol=1.e-12)                                                     
      common/comsic/ dt,gvdt(km),svdt(km),amdt(km*km),bmdt(km*km),              
     1               dm(km*km,0:jcap)                                           
#define DEFAULT
#ifdef MINV
#undef DEFAULT
      dimension  work(2*km)                                                     
#endif
#ifdef DEFAULT
      dimension iwork(2*km)                                                     
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      dt=deltim                                                                 
      do k=1,km                                                                 
        gvdt(k)=dt*gv(k)                                                        
        svdt(k)=dt*sv(k)                                                        
      enddo                                                                     
      do kj=1,km*km                                                             
        amdt(kj)=dt*am(kj)                                                      
        bmdt(kj)=dt*bm(kj)                                                      
        dm(kj,0)=0.                                                             
      enddo                                                                     
      do kj=1,km*km,km+1                                                        
        dm(kj,0)=1.                                                             
      enddo                                                                     
#ifdef ORIGIN_THREAD
c$doacross share(dt,dm,cm,km,tol),
c$&        local(n,dt2nn1,kj,work,iwork,det)                                    
#endif
#ifdef CRAY_THREAD
cmic$ do all shared(dt,dm,cm,km,tol) private(n,dt2nn1,kj,work,iwork,det)
#endif
#ifdef OPENMP
!$omp parallel do private(n,dt2nn1,kj,work,iwork,det)
#endif
      do n=1,jcap                                                               
        dt2nn1=dt**2*(n*(n+1))                                                  
        do kj=1,km*km                                                           
#ifndef HYBRID
          dm(kj,n)=dm(kj,0)-dt2nn1*cm(kj)
#else
          dm(kj,n)=dm(kj,0)+dt2nn1*cm(kj)
#endif
        enddo                                                                   
#define DEFAULT
#ifdef MINV
#undef DEFAULT
        call  minv(dm(1,n),km,km,work,det,tol,0,1)                          
#endif
#ifdef DEFAULT
        call iminv(dm(1,n),km,det,iwork(1),iwork(km+1))                       
#endif
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
