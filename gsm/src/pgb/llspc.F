#define IO io
#define JO jo
#define KO ko
#include <define.h>
      subroutine llspc(im2,jm2,km,m,ncpus,clat,slat,wlat,trig,ifax,             
     &                 u,v,o,z,t,a)                                             
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram: llspc          perform spectral functions on latlon fields        
c   prgmmr: iredell          org: w/nmc23     date: 94-08-31                    
c                                                                               
c abstract: this subprogram performs spectral functions on latlon fields        
c   winds, vertical velocities, heights and temperatures are smoothed           
c   and the absolute vorticity is computed from the winds.                      
c   all fields are on a global latlon grid that includes the poles              
c   with northern and southern hemispheres paired.                              
c                                                                               
c program history log:                                                          
c   95-08-31  iredell                                                           
c                                                                               
c usage:    call llspc(im2,jm2,km,m,ncpus,clat,slat,wlat,trig,ifax,             
c    &                 u,v,o,z,t,a)                                             
c   input arguments:                                                            
c     im2          integer twice the number of longitude points                 
c     jm2          integer half the number of latitude points                   
c     km           integer number of levels                                     
c     m            integer spectral truncation (jm2-2 recommended)              
c     ncpus        integer number of cpus over which to distribute work         
c     clat         real (jm2) cosines of latitude                               
c     slat         real (jm2) positive sines of latitude                        
c     wlat         real (jm2) gaussian weights                                  
c     trig         real (im2) trigonometric quantities for the fft              
c     ifax         integer (20) factors for the fft                             
c     u            real (im2,jm2,km) zonal wind in m/s                          
c     v            real (im2,jm2,km) meridional wind in m/s                     
c     o            real (im2,jm2,km) vertical velocity in pa/s                  
c     z            real (im2,jm2,km) heights in m                               
c     t            real (im2,jm2,km) temperature in k                           
c   output arguments:                                                           
c     u            real (im2,jm2,km) smoothed zonal wind in m/s                 
c     v            real (im2,jm2,km) smoothed meridional wind in m/s            
c     o            real (im2,jm2,km) smoothed vertical velocity in pa/s         
c     z            real (im2,jm2,km) smoothed heights in m                      
c     t            real (im2,jm2,km) smoothed temperature in k                  
c     a            real (im2,jm2,km) absolute vorticity in 1/s                  
c                                                                               
c subprograms called:                                                           
c   gspc         compute spectral constants                                     
c   rfftmlt      compute fft                                                    
c   pleg         compute legendre functions                                     
c   panaly       compute spectral from fourier                                  
c   uv2dz        compute divergence and vorticity in spectral space             
c   dz2uv        compute wind components in spectral space                      
c   psynth       compute fourier from spectral                                  
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      real slat(jm2),clat(jm2),wlat(jm2)                                        
      real trig(im2)                                                            
      integer ifax(20)                                                          
      real u(im2,jm2,km),v(im2,jm2,km),o(im2,jm2,km)                            
      real z(im2,jm2,km),t(im2,jm2,km),a(im2,jm2,km)                            
c                                                                               
#define DEFAULT
#ifdef DYNAMIC_ALLOC
#undef DEFAULT
      real eps((m+1)*(m+2)/2),epstop(m+1)                                       
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
      real eon((m+1)*(m+2)/2),eontop(m+1)                                       
      real wfft(im2,2*6*km)                                                     
      real pln((m+1)*(m+2)/2),plntop(m+1)                                       
      real f1(im2/2+3,2,5*km,ncpus),f2(im2/2+3,2,6*km)                          
      real s1((m+1)*(m+2)+1,6*km),s1top(2*(m+1),6*km)                           
      real sd((m+1)*(m+2))                                                      
      integer mp(6*km)                                                          
#endif
      parameter(mfil=(jo_+1)/2-2)                                               
#ifdef DEFAULT
      real eps((mfil+1)*(mfil+2)/2),epstop(mfil+1)                              
      real enn1((mfil+1)*(mfil+2)/2),elonn1((mfil+1)*(mfil+2)/2)                
      real eon((mfil+1)*(mfil+2)/2),eontop(mfil+1)                              
      real wfft(2*io_,2*6*ko_)                                                  
      real pln((mfil+1)*(mfil+2)/2),plntop(mfil+1)                              
      parameter(mcpus=1)                                                        
      real f1(2*io_/2+3,2,5*ko_,mcpus),f2(2*io_/2+3,2,6*ko_)                    
      real s1((mfil+1)*(mfil+2)+1,6*ko_),s1top(2*(mfil+1),6*ko_)                
      real sd((mfil+1)*(mfil+2))                                                
      integer mp(6*ko_)                                                         
#endif
c                                                                               
      parameter(omega= 7.2921e-5 )                                              
c                                                                               
c  set transform constants                                                      
c                                                                               
      im=im2/2                                                                  
      ix=im+3                                                                   
c                                                                               
      do k=1,km*2                                                               
        mp(k)=1                                                                 
      enddo                                                                     
      do k=2*km+1,6*km                                                          
        mp(k)=0                                                                 
      enddo                                                                     
      if(m.ne.mfil) then                                                        
        write(6,*) 'm.ne.nfil in llspc'                                         
        call abort                                                              
      endif                                                                     
c                                                                               
      nc=(m+1)*(m+2)+1                                                          
      nctop=2*(m+1)                                                             
c                                                                               
      call gspc(m,eps,epstop,enn1,elonn1,eon,eontop)                            
c                                                                               
c  transform to spectral space                                                  
c                                                                               
      do k=1,km                                                                 
        do i=1,nc                                                               
          s1(i,k)=0.                                                            
          s1(i,k+km)=0.                                                         
          s1(i,k+2*km)=0.                                                       
          s1(i,k+3*km)=0.                                                       
          s1(i,k+4*km)=0.                                                       
        enddo                                                                   
        do i=1,nctop                                                            
          s1top(i,k)=0.                                                         
          s1top(i,k+km)=0.                                                      
          s1top(i,k+2*km)=0.                                                    
          s1top(i,k+3*km)=0.                                                    
          s1top(i,k+4*km)=0.                                                    
          s1top(i,k+5*km)=0.                                                    
        enddo                                                                   
      enddo                                                                     
c                                                                               
c     warning!! wrong in original code 'do j1=2,jm2,ncpus'                      
c                                                                               
      do j1=1,jm2,ncpus                                                         
        j2=min(j1+ncpus-1,jm2)                                                  
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 shared(j1,j2,km,im,f1,clat,u,v,o,z,t,trig,ifax,ix)                       
cmic$1 private(j,jc,i,k,wfft)                                                   
#endif
        do j=j1,j2                                                              
          jc=j-j1+1                                                             
          do k=1,km                                                             
            do i=1,im                                                           
              if(clat(j).ne.0.) then                                            
                f1(i,1,k,jc)=u(i,j,k)/clat(j)**2                                
                f1(i,2,k,jc)=u(i+im,j,k)/clat(j)**2                             
                f1(i,1,k+km,jc)=v(i,j,k)/clat(j)**2                             
                f1(i,2,k+km,jc)=v(i+im,j,k)/clat(j)**2                          
              else                                                              
                f1(i,1,k,jc)=0.                                                 
                f1(i,2,k,jc)=0.                                                 
                f1(i,1,k+km,jc)=0.                                              
                f1(i,2,k+km,jc)=0.                                              
              endif                                                             
              f1(i,1,k+2*km,jc)=o(i,j,k)                                        
              f1(i,2,k+2*km,jc)=o(i+im,j,k)                                     
              f1(i,1,k+3*km,jc)=z(i,j,k)                                        
              f1(i,2,k+3*km,jc)=z(i+im,j,k)                                     
              f1(i,1,k+4*km,jc)=t(i,j,k)                                        
              f1(i,2,k+4*km,jc)=t(i+im,j,k)                                     
            enddo                                                               
          enddo                                                                 
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
            call rfftmlt(f1(1,1,1,jc),wfft,trig,ifax,1,ix,im,2*5*km,-1)         
#endif
#ifdef DEFAULT
          call fft99m (f1(1,1,1,jc),wfft,trig,ifax,1,ix,im,2*5*km,-1)           
#endif
        enddo                                                                   
        do j=j1,j2                                                              
          jc=j-j1+1                                                             
          call pleg(m,slat(j),clat(j),eps,epstop,pln,plntop)                    
          wa=wlat(j)                                                            
          call panaly(m,im,ix,nc,nctop,5*km,wa,clat(j),pln,plntop,mp,           
     &                f1(1,1,1,jc),s1,s1top)                                    
        enddo                                                                   
      enddo                                                                     
c                                                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  spectrally compute vorticity                                                 
c                                                                               
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 shared(km,m,enn1,elonn1,eon,eontop,s1,s1top)                             
cmic$1 private(k,sd)                                                            
#endif
c                                                                               
      do k=1,km                                                                 
        call uv2dz(m,enn1,elonn1,eon,eontop,                                    
     &             s1(1,k),s1(1,k+km),s1top(1,k),s1top(1,k+km),                 
     &             sd,s1(1,k+5*km))                                             
        call dz2uv(m,enn1,elonn1,eon,eontop,                                    
     &             sd,s1(1,k+5*km),                                             
     &             s1(1,k),s1(1,k+km),s1top(1,k),s1top(1,k+km))                 
      enddo                                                                     
c                                                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  transform to grid                                                            
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 shared(jm2,m,slat,clat,eps,epstop)                                       
cmic$1 shared(im,ix,nc,nctop,km,mp,s1,s1top,trig,ifax,u,v,o,z,t,a)              
cmic$1 private(j,pln,plntop,f2,wfft,k,i)                                        
#endif
c                                                                               
      do j=1,jm2                                                                
        call pleg(m,slat(j),clat(j),eps,epstop,pln,plntop)                      
        call psynthv(m,im,ix,nc,nctop,6*km,clat(j),pln,plntop,mp,               
     &              s1,s1top,f2)                                                
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
          call rfftmlt(f2,wfft,trig,ifax,1,ix,im,2*6*km,1)                      
#endif
#ifdef DEFAULT
        call fft99m (f2,wfft,trig,ifax,1,ix,im,2*6*km,1)                        
#endif
        do k=1,km                                                               
          do i=1,im                                                             
            u(i,j,k)=f2(i,1,k)                                                  
            u(i+im,j,k)=f2(i,2,k)                                               
            v(i,j,k)=f2(i,1,k+km)                                               
            v(i+im,j,k)=f2(i,2,k+km)                                            
            o(i,j,k)=f2(i,1,k+2*km)                                             
            o(i+im,j,k)=f2(i,2,k+2*km)                                          
            z(i,j,k)=f2(i,1,k+3*km)                                             
            z(i+im,j,k)=f2(i,2,k+3*km)                                          
            t(i,j,k)=f2(i,1,k+4*km)                                             
            t(i+im,j,k)=f2(i,2,k+4*km)                                          
            a(i,j,k)=f2(i,1,k+5*km)+2*omega*slat(j)                             
            a(i+im,j,k)=f2(i,2,k+5*km)-2*omega*slat(j)                          
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
