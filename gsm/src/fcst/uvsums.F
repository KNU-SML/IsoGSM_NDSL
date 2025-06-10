#define JCAP1 jcap1
#define JCAP jcap
#define JCAPP jcapp

      subroutine uvsums(fpu,fmu,fpv,fmv,topuln,topvln,qvv,wgt,
     1                  llstr,llens,lwvdef,levs)             
c                                                                               
c            this sr assumes jcap is even                                       
c            this sr assumes jcap is even                                       
c            this sr assumes jcap is even                                       
c                                                                               
#include <paramodel.h>
#ifdef MP
#define JCAPS jcapp_
#else
#define JCAPS jcap_
#endif
      dimension     fpu(2,0:JCAPS,levs)
      dimension     fpv(2,0:JCAPS,levs)
      dimension     fmu(2,0:JCAPS,levs)
      dimension     fmv(2,0:JCAPS,levs)
      dimension   topuln(2,0:JCAPS,levs)
      dimension   topvln(2,0:JCAPS,levs)
      dimension   toppln(0:JCAPS)                                               
      dimension   qvv(2,0:JCAPS)                                                   
      dimension   lwvdef(jcap1_)
c                                                                               
c     ----------------------------------------------------------------          
c     compute expansion coeffs. for top rows of u and v                         
c     ----------------------------------------------------------------          
cc                                                                              
cc                                                                              
      do 10 l=0,JCAPS                                                           
      toppln(l) = qvv(1,l)*wgt                                                    
   10 continue                                                                  
cc                                                                              
cc
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 autoscope                                                                
#endif
c    
      do lx=1,llens
        l=lwvdef(llstr+lx)
        ll=lx-1
c
        if( mod(l,2).eq.1 ) then
c                                                                               
c  even wavenumber for l=1,3,5,7,.........
c
c         compute the even (n-l) expansion coefficients for each level           
c         ------------------------------------------------------------           
c                                                                               
c         real part                                                              
c                                                                               
          do k=1,levs
            topuln(1,ll,k) = topuln(1,ll,k)+fpu(1,ll,k)*toppln(ll)
            topvln(1,ll,k) = topvln(1,ll,k)+fpv(1,ll,k)*toppln(ll)
          enddo
c                                                                               
c         imaginary part                                                         
c                                                                               
          do k=1,levs
            topuln(2,ll,k) = topuln(2,ll,k)+fpu(2,ll,k)*toppln(ll)
            topvln(2,ll,k) = topvln(2,ll,k)+fpv(2,ll,k)*toppln(ll)
          enddo
c
c                                                                               
c  odd wavenumber for l=0,2,4,6,.........
c
        else
c
c         compute the odd (n-l) expansion coefficients for each level            
c         -----------------------------------------------------------            
c
c                                                                               
c         real part                                                              
c
          do k=1,levs
            topuln(1,ll,k) = topuln(1,ll,k)+fmu(1,ll,k)*toppln(ll)
            topvln(1,ll,k) = topvln(1,ll,k)+fmv(1,ll,k)*toppln(ll)
          enddo
c                                                                               
c         imaginary part                                                         
c                                                                               
          do k=1,levs
            topuln(2,ll,k) = topuln(2,ll,k)+fmu(2,ll,k)*toppln(ll)
            topvln(2,ll,k) = topvln(2,ll,k)+fmv(2,ll,k)*toppln(ll)
          enddo
c
c
        endif
c
      enddo                                                                    
c
      return                                                                    
      end                                                                       
