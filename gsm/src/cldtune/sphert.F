#define JCAP jcap
#define LATG latg
#define LONF lonf
#include <define.h>
      subroutine sphert(idir,grid,wave,mlt,fac,imax,jmax,maxwv,iromb)           
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
