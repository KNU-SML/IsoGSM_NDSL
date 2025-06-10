#define IO io
#define JCAP jcap
#define JO jo
#define LEVS levs
#define RD rd
#define RV rv
#include <define.h>
#ifdef SCM
      subroutine trss(nflds,ntrace,ncloud,ss,f)
#else
      subroutine trss(nflds,ntrace,ncloud,trig,ifax,eps,epstop,ss,sstop,
     &                coslat,sinlat,f)                                          
#endif
#include <paramodel.h>
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)                         
      parameter(jcap=jcap_,levs=levs_)                                          
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)                     
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    trss        transform spectral to grid                         
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: transforms spectral to gridded data on a latitude pair              
c           and computes dry temperature and surface pressure                   
c           and winds and gradients without a cosine latitude factor.           
c           subprogram rdss should be called already                            
c           to read spectral data and initialize utility fields.                
c           this subprogram can be called from a multiprocessed segment.        
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call trss(jcap,nc,nctop,levs,trig,ifax,eps,epstop,ss,sstop,         
c    &                io2,io22,coslat,sinlat,f)                                 
c                                                                               
c   input argument list:                                                        
c     jcap     - integer spectral truncation                                    
c     nc       - integer number of spectral coefficients                        
c     nctop    - integer number of spectral coefficients over top               
c     levs     - integer number of levels                                       
c     trig     - real (io2) trigonometric quantities for the fft                
c     ifax     - integer (20) factors for the fft                               
c     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))        
c     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top            
c     ss       - real (nc,6*levs+6) spectral coefs                              
c     sstop    - real (nctop,6*levs+6) spectral coefs over top                  
c     io2    - integer number of valid data points per latitude pair            
c     io22   - integer longitude dimension of data (>=io2+4)                    
c     coslat   - real cosine of latitude of the latitude pair                   
c     sinlat   - real sine of latitude of the northern latitude                 
c                                                                               
c   output argument list:                                                       
c     f        - real (io22,6*levs+6) gridded data                              
c                (:,1:levs)             vorticity                               
c                (:,levs+1:2*levs)      divergence                              
c                (:,2*levs+1:3*levs)    temperature                             
c                (:,3*levs+1:4*levs)    specific humidity                       
c                (:,4*levs+1)           d(lnps)/dx                              
c                (:,4*levs+2)           d(lnps)/dy                              
c                (:,4*levs+3:5*levs+2)  zonal wind                              
c                (:,5*levs+3:6*levs+2)  meridional wind                         
c                (:,6*levs+3)           surface pressure                        
c                (:,6*levs+4)           orography                               
c                (:,6*levs+5)           d(orog)/dx                              
c                (:,6*levs+6)           d(orog)/dy                              
c                                                                               
c subprograms called:                                                           
c   pleg         compute associated legendre polynomials                        
c   psynth       synthesize fourier from spectral coefficients                  
c   rfftmlt      fast fourier transform                                         
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#ifdef SCM
      real ss(nc,nflds),f(io22,nflds)
#else
      dimension trig(io2),ifax(20)                                              
      real eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)                              
      real ss(nc,nflds),sstop(nctop,nflds),f(io22,nflds)
      real pln((jcap+1)*(jcap+2)/2),plntop(jcap+1)                              
      real wfft(io22,2*(5*levs+ntrace*levs+6+ncloud*levs))
      integer mp(5*levs+ntrace*levs+6+ncloud*levs)
#endif
      parameter(fv= rv_ / rd_ -1.)                                  
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef SCM
      print*,' trss ',nflds,ntrace,ncloud
      do n=1,nflds
        do i=1,io22
          f(i,n)=0.0    ! initialization
        enddo
      enddo
      do n=1,nflds
        do i=1,io2
          f(i,n)=ss(i,n)
        enddo
      enddo
#else
c  transform spectral coefficients to fourier coefficients                      
      print*,' trss ',nflds,ntrace,ncloud,trig,ifax,coslat,sinlat
      call pleg(jcap,sinlat,coslat,eps,epstop,pln,plntop)                       
      call psynth(jcap,io22/2,nc,nctop,nflds,pln,plntop,ss,sstop,f)          
      print*,' trss ',(f(1,k),k=nflds-levs_+1,nflds)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  transform fourier coefficients to gridded data                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(f,wfft,trig,ifax,1,io22/2,io2/2,2*nflds,1)            
#endif
#ifdef DEFAULT
      call fft99m (f,wfft,trig,ifax,1,io22/2,io2/2,2*nflds,1)              
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  move southern hemisphere latitude after northern hemisphere latitude         
      do k=1,nflds
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
        do i=1,io2/2                                                            
          f(io2/2+i,k)=f(io22/2+i,k)                                            
        enddo                                                                   
      enddo                                                                     
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute dry temperature from virtual temperature                             
c  and surface pressure from log surface pressure                               
c  and divide gradients and winds by cosine of latitude.                        
      do k=1,levs                                                               
        do i=1,io2                                                              
          f(i,4*levs+k+6)=f(i,4*levs+k+6)/(1.+fv*f(i,5*levs+k+6))                     
#ifndef SCM
          if(coslat.ne.0.) then                                                 
            f(i,2*levs+6+k)=f(i,2*levs+6+k)/coslat                              
            f(i,3*levs+6+k)=f(i,3*levs+6+k)/coslat                              
          else                                                                  
            f(i,2*levs+6+k)=0.                                                  
            f(i,3*levs+6+k)=0.                                                  
          endif                                                                 
#endif
        enddo                                                                   
      enddo                                                                     
      do i=1,io2                                                                
        f(i,4)=exp(f(i,4))                                        
      enddo                                                                     
#ifndef SCM
      if(coslat.ne.0.) then                                                     
        do i=1,io2                                                              
          f(i,5)=f(i,5)/coslat                                    
          f(i,6)=f(i,6)/coslat                                    
          f(i,2)=f(i,2)/coslat                                    
          f(i,3)=f(i,3)/coslat                                    
        enddo                                                                   
      else                                                                      
        do i=1,io2                                                              
          f(i,5)=0.                                                      
          f(i,6)=0.                                                      
          f(i,2)=0.                                                      
          f(i,3)=0.                                                      
        enddo                                                                   
      endif                                                                     
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      print*,' trss ',(f(1,k),k=nflds-levs_+1,nflds)
      return                                                                    
      end                                                                       
