#define JCAP jcap
#define LATG latg
#define LONF lonf
#include <define.h>
      subroutine sphert(idir,grid,wave,mlt,fac,imax,jmax,maxwv,iromb)           
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:  sphert     spherical transform                                   
c   prgmmr: kanamitsu     org: w/nmc23       date: 92-04-16                     
c                                                                               
c abstract: transforms a field between grid and spectral domains.               
c   this versatile routine will:                                                
c     ...transform grid to spectral or transform spectral to grid;              
c     ...pass a gaussian grid or an equally-space grid;                         
c     ...pass a triangular truncation or a rhomboidal truncation;               
c     ...optionally transform with derivatives of legendre functions;           
c     ...optionally transform from grid dividing by coslat**2;                  
c     ...optionally multiply spectral field by complex factors.                 
c                                                                               
c program history log:                                                          
c   92-04-16  iredell                                                           
c                                                                               
c usage:    call sphert(idir,grid,wave,mlt,fac,imax,jmax,maxwv,iromb)           
c                                                                               
c   input argument list:                                                        
c     idir     - must be one of 1,2,3,4,101,102,103,104,-1,-2,-101,-102,        
c                idir> 0 to transform grid to wave,                            
c                idir < 0 to transform wave to grid,                            
c                abs(idir) < 100 for gaussian grid,                             
c                abs(idir)> 100 for equally-spaced grid,                       
c                idir is odd for normal transform,                              
c                idir is even for legendre derivative transform,                
c                idir last digit> 2 to transform with cosine latitude.         
c     grid     - if idir> 0, real (imax,jmax) field to transform.              
c                grid starts at north pole and greenwich meridian.              
c     wave     - if idir < 0, complex (kmax) field to transform,                
c                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
c                wave starts at the global mean component and then              
c                contains the zonally symmetric components.                     
c     mlt      - multiplication option on wave values.                          
c                mlt = 0 for no multiplication,                                 
c                mlt = 1 to multiply wave by fac                                
c                mlt = -1 to multiply wave by sqrt(-1)*fac                      
c     fac      - if mlt.ne.0, real (kmax) field to multiply wave.               
c     imax     - longitudinal dimension of the grid                             
c     jmax     - latitudinal dimension of the grid                              
c     maxwv    - spectral truncation of the wave                                
c     iromb    - iromb = 0 for triangular truncation                            
c                iromb = 1 for rhomboidal truncation                            
c                                                                               
c   output argument list:                                                       
c     grid     - if idir < 0, real (imax,jmax) field output.                    
c                grid starts at north pole and greenwich meridian.              
c     wave     - if idir> 0, complex (kmax) field output,                      
c                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
c                wave starts at the global mean component and then              
c                contains the zonally symmetric components.                     
c                                                                               
c   subprograms called:                                                         
c     unique:                                                                   
c     gausslat   - compute gaussian latitudes                                   
c     equallat   - compute equally-spaced latitudes                             
c     (legendre) - compute legendre polynomials                                 
c     fftfax     - fft (library call can be substituted)                        
c     rfftmlt    - fft (library call can be substituted)                        
c                                                                               
c   remarks: fortran 9x extensions are used.                                    
c     fpp can be used to inline subprograms.                                    
c     transforming with an equally-spaced grid is irreversible.                 
c                                                                               
c attributes:                                                                   
c   cray ymp.                                                                   
c                                                                               
c$$$                                                                            
#ifdef CRAY_THREAD
cfpp$ expand(legendre)                                                          
#endif
#include <paramodel.h>
      save                                                                      
      real grid(imax,jmax)                                                      
      complex wave((maxwv+1)*((iromb+1)*maxwv+2)/2)                             
      real fac((maxwv+1)*((iromb+1)*maxwv+2)/2)                                 
c                                                                               
      parameter(ngg=2*(1+(lonf_+1)/2))                                          
      real gg(ngg,latg_)                                                        
      parameter(npnm=(jcap_+1)*((0+1)*jcap_+2)/2)                               
      complex ww(npnm),w2(-1:1),ws                                              
      parameter(ntrigs=2*lonf_)                                                 
c-we  real trigs(ntrigs),ifax(20)                                               
      real trigs(ntrigs)                                                        
      integer ifax(20)                                                          
c                                                                               
      real cosclt(latg_),wgtclt(latg_)                                          
      parameter(nwrkft=2*lonf_*latg_)                                           
      real wrkfft(nwrkft)                                                       
      integer is(npnm)                                                          
      real pnm  (npnm)                                                          
      parameter(nepx=jcap_+1)                                                   
      real ex(0:nepx),px(-1:nepx)                                               
c                                                                               
      kmax=(maxwv+1)*((iromb+1)*maxwv+2)/2                                      
      ipd=1-mod(abs(idir),2)                                                    
      isd=1-2*ipd                                                               
      icd=(mod(abs(idir),10)-1)/2                                               
      jump=2*(1+(imax+1)/2)                                                     
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
            call fftfax (lonf_,ifax,trigs)                                      
#endif
#ifdef DEFAULT
            call    fax (ifax, lonf_,3)                                         
            call fftrig (trigs,lonf_,3)                                         
#endif
      if(abs(idir).lt.100) then                                                 
        call gausslat(jmax,cosclt,wgtclt)                                       
      else                                                                      
        call equallat(jmax,cosclt,wgtclt)                                       
      endif                                                                     
      k=0                                                                       
      do m=0,maxwv                                                              
        do n=m,iromb*m+maxwv                                                    
          k=k+1                                                                 
          is(k)=isd*(1-2*mod(n-m,2))                                            
        enddo                                                                   
      enddo                                                                     
                                                                                
      if(idir.gt.0) then                                                        
                                                                                
        do j=1,jmax                                                             
          do i=1,imax                                                           
            gg(i,j)=grid(i,j)                                                   
          enddo                                                                 
        enddo                                                                   
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)                  
#endif
#ifdef DEFAULT
        call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)                  
#endif
        do k=1,kmax                                                             
          wave(k)=0.                                                            
        enddo                                                                   
        do j=1,(jmax+1)/2                                                       
          jr=jmax+1-j                                                           
          call legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)                    
          wj=wgtclt(j)                                                          
          if(icd.ne.0.and.cosclt(j).lt.1.) wj=wj/(1.-cosclt(j)**2)              
          k=0                                                                   
          do m=0,maxwv                                                          
            w2(1)=wj*cmplx(gg(2*m+1,j)+gg(2*m+1,jr),                            
     &                     gg(2*m+2,j)+gg(2*m+2,jr))                            
            w2(-1)=wj*cmplx(gg(2*m+1,j)-gg(2*m+1,jr),                           
     &                      gg(2*m+2,j)-gg(2*m+2,jr))                           
            do n=m,iromb*m+maxwv                                                
              k=k+1                                                             
              wave(k)=wave(k)+w2(is(k))*pnm(k)                                  
            enddo                                                               
          enddo                                                                 
        enddo                                                                   
        if(mlt.ne.0) then                                                       
          ws=cmplx(1.,0.)                                                       
          if(mlt.lt.0) ws=cmplx(0.,1.)                                          
          do k=1,kmax                                                           
            wave(k)=wave(k)*ws*fac(k)                                           
          enddo                                                                 
        endif                                                                   
                                                                                
      else                                                                      
                                                                                
        if(mlt.ne.0) then                                                       
          ws=cmplx(1.,0.)                                                       
          if(mlt.lt.0) ws=cmplx(0.,1.)                                          
          do k=1,kmax                                                           
            ww(k)=wave(k)*ws*fac(k)                                             
          enddo                                                                 
        else                                                                    
          do k=1,kmax                                                           
            ww(k)=wave(k)                                                       
          enddo                                                                 
        endif                                                                   
        do j=1,(jmax+1)/2                                                       
          call legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)                    
          k=0                                                                   
          do m=0,maxwv                                                          
            w2(1)=0.                                                            
            w2(-1)=0.                                                           
            do n=m,iromb*m+maxwv                                                
              k=k+1                                                             
              w2(1)=w2(1)+ww(k)*pnm(k)                                          
              w2(-1)=w2(-1)+ww(k)*pnm(k)*is(k)                                  
            enddo                                                               
            gg(2*m+1,jmax+1-j)=real(w2(-1))                                     
            gg(2*m+2,jmax+1-j)=aimag(w2(-1))                                    
            gg(2*m+1,j)=real(w2(1))                                             
            gg(2*m+2,j)=aimag(w2(1))                                            
          enddo                                                                 
        enddo                                                                   
        do i=2*maxwv+3,jump                                                     
          do j=1,jmax                                                           
            gg(i,j)=0.                                                          
          enddo                                                                 
        enddo                                                                   
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)                   
#endif
#ifdef DEFAULT
        call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)                   
#endif
        do j=1,jmax                                                             
          do i=1,imax                                                           
            grid(i,j)=gg(i,j)                                                   
          enddo                                                                 
        enddo                                                                   
                                                                                
      endif                                                                     
                                                                                
      return                                                                    
      end                                                                       
