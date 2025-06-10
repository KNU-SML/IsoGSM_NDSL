#define G g
#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define RD rd
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine newps(t,ps,hold,hnew,psnew,sip)
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    newps       interpolate surface pressure.                      
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: using the hydrostatic equation, the surface pressure                
c           is interpolated from input orography to output orography.           
c           below the input surface, the temperature lapse rate                 
c           is fixed at -6.5k/km.  above the input surface,                     
c           surface pressure is interpolated from height profiles               
c           computed as in the mrf model.                                       
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell                                                      
c                                                                               
c usage:    call newps(t,q0,h0,h1,q1,si,sl)                                     
c   input argument list:                                                        
c     t        - temperature                                                    
c     q0       - old ln(psfc)                                                   
c     h0       - old orography                                                  
c     h1       - new orography                                                  
c     si       - edge sigma values                                              
c     sl       - full sigma values                                              
c                                                                               
c   output argument list:                                                       
c     q1       - new ln(psfc)                                                   
c                                                                               
c   subprograms called:                                                         
c     bnmc     - compute mrf temperature to height matrix                       
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
c                                                                               
#include <paramter.h>
#include <parmchgr.h>
c                                                                               
      parameter(idimt=idim*2,jdimhf=jdim/2,                                     
     6          kdimp=kdim+1,kdimm=kdim-1,                                      
     7          kdimt=kdim*2)                                                   
c                                                                               
         dimension t(idimt,kdimi),ps(idimt),psnew(idimt)                        
         dimension hold(idimt),hnew(idimt),sip(idimt,kdimi+1)
         dimension tau(idimt,kdimi)                                             
         dimension rlsig(idimt,kdimi+1)
         dimension h(idimt,kdimi+1)                                             
         dimension a(idimt)                                                     
         dimension absa(idimt)                                                  
c--------                                                                       
c-------- compute interface heights.                                            
c--------                                                                       
c        print *,' ------ into newpsfc ---'                                     
         rog = rd_/g_                                                           
         do k=1,kdimi                                                         
           do i=1,idimt
             rlsig(i,k) = - alog(sip(i,k))
             tau(i,k) = t(i,k) * rog
           enddo
         enddo 
         do 6 i=1,idimt
         rlsig(i,kdimi+1) = - alog ( .005 )
6        continue
c.........                                                                      
         do 10 i=1,idimt                                                        
         psnew(i)= - ps(i)                                                      
         h(i,1)=hold(i)                                                         
10       continue                                                               
c.........                                                                      
         do 100 k=1,kdimi                                                       
         do 100 i=1,idimt                                                       
           h(i,k+1)=h(i,k)+tau(i,k)*(rlsig(i,k+1)-rlsig(i,k))
100      continue                                                               
c--------                                                                       
c-------- loop over layers, testing to see if new sfc pressure is               
c-------- in layer, and obtaining it, if so.                                    
c--------                                                                       
         eps=1.e-1                                                             
         do 200 k=1,kdimi                                                       
c----------                                                                     
c---------- compute lapse rate                                                  
c----------                                                                     
           kp=min(k+1,kdimi)                                                    
           km=max(1,k-1)                                                        
           kppz=min(k+2,kdimi)                                                  
           kpz=kppz-1                                                           
           kmmz=max(1,k-2)                                                      
           kmz=kmmz+1                                                           
           do 2020 iq=1,idimt                                                   
           a(iq)=2.0e0*(tau(iq,kp)-tau(iq,km))/                                
     1          (rlsig(iq,kpz)+rlsig(iq,kppz)
     2          -rlsig(iq,kmz)-rlsig(iq,kmmz))
           absa(iq)=abs(a(iq))                                                  
2020       continue                                                             
           do 2030 i=1,idimt                                                    
           if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                    
     1        absa(i).gt.eps) psnew(i)=rlsig(i,k)
     1                                -ps(i)+(sqrt(tau(i,k)**2+                 
     2        2.0e0*a(i)*(hnew(i)-h(i,k))) - tau(i,k) )/a(i)                   
c.........                                                                      
           if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                    
     1        absa(i).le.eps) psnew(i)=rlsig(i,k)
     2                                -ps(i)+(hnew(i)-h(i,k))/tau(i,k)          
2030       continue                                                             
200      continue                                                               
c--------                                                                       
c-------- do points which fall below first layer (use fixed lapse rate          
c--------  of 6.5 deg per km.                                                   
c--------                                                                       
         gamma=6.5e-3                                                           
         gascon=rd_                                                             
         g=g_                                                                   
         c=gascon*gamma/g                                                       
         do 8010 iq=1,idimt                                                     
         a(iq)=tau(iq,1)*
     1       (1.e0-exp(c*(rlsig(iq,1)-rlsig(iq,2))))
     2                  /(rlsig(iq,1)-rlsig(iq,2)) 
8010     continue                                                               
         do 8050 i=1,idimt                                                      
         if( hnew(i).le.h(i,1) ) psnew(i)=rlsig(i,1)
     1                                   -ps(i)+(sqrt(tau(i,1)**2+              
     1     2.e0*a(i)*(hnew(i)-h(i,1))) - tau(i,1))/a(i)                        
8050     continue                                                               
         icount = 0                                                             
         do 8100 i=1,idimt                                                      
         psnew(i) = - psnew(i)                                                  
ccc      if( psnew(i).eq.ps(i) ) icount = icount + 1                            
         if( psnew(i).eq.ps(i) ) then                                           
              icount = icount + 1                                               
         else                                                                   
c         print *,'  gz gzold ps psnew ',hold(i),hnew(i),ps(i),psnew(i)         
         endif                                                                  
8100     continue                                                               
c                                                                               
c        write(6,300)icount                                                     
300      format(' number of surface pressure points not updated=',i6)           
c                                                                               
       return                                                                   
       end                                                                      
