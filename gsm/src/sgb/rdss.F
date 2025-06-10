#define IOSGB iosgb
#define JCAP jcap
#define JOSGB josgb
#define LEVS levs
#include <define.h>
      subroutine rdss(nss,si,sl,idrt,                                              
     $                mflds,ntrace,ncloud,
     1                clat,slat,trig,ifax,eps,epstop,ss,sstop)                  
c                                                                               
#include <paramodel.h>
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)                         
      parameter(jcap=jcap_,levs=levs_)                                          
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)                     
c
c  number of variables for spectral conversion
c
c      parameter(mflds=13*levs_+6)                                                
      dimension si(levs+1),sl(levs),clat(johf),slat(johf),
     1          trig(io2),ifax(20)               
      real eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)                              
      real ss(nc,mflds),sstop(nctop,mflds)                                      
c                                                                               
      real wlat(johf)                                                           
      real enn1((jcap+1)*(jcap+2)),elonn1((jcap+1)*(jcap+2)/2)                  
      real eon((jcap+1)*(jcap+2)/2),eontop(jcap+1)                              
c                                                                               
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 )                                   
c                                                                               
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    rdss        read data from a sigma spectral file               
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: reads the records of orography, surface pressure,                   
c           divergence and vorticity, temperature and humidity                  
c           from a sigma spectral file.  it is assumed that the first           
c           two header records of the file have already been read.              
c           the gradients of orography and log surface pressure                 
c           and the wind components are also computed in spectral space.        
c           the geopotential of the pressure gradient is computed too.          
c           also, some spectral transform utility fields are computed.          
c           subprogram trss should be used to transform to grid                 
c           as well as compute dry temperature and surface pressure             
c           and winds and gradients without a cosine latitude factor.           
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call rdss(nss,jcap,nc,nctop,johf,io2,levs,sl,                       
c   &                 clat,slat,trig,ifax,eps,epstop,ss,sstop)                  
c                                                                               
c   input argument list:                                                        
c     nss      - integer unit from which to read file                           
c     sl       - real (levs) sigma full level values                            
c     si       - real (levs) sigma interface level values                            
c                                                                               
c   output argument list:                                                       
c     clat     - real (johf) cosines of latitude                                
c     slat     - real (johf) sines of latitude                                  
c     trig     - real (io2) trigonometric quantities for the fft                
c     ifax     - integer (20) factors for the fft                               
c     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))        
c     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top            
c     ss       - real (nc,6*levs+6) spectral coefs                              
c     sstop    - real (nctop,6*levs+6) spectral coefs over top                  
c
c            (:, 1:levs              ) vorticity
c            (:, 1*levs+ 1: 2*levs   ) divergence
c            (:, 2*levs+ 1: 3*levs   ) temperature
c            (:, 3*levs+ 1: 4*levs   ) specific humidity
c            (:, 4*levs+ 1: 5*levs   ) u-component of wind
c            (:, 5*levs+ 1: 6*levs   ) v-component of wind
c            (:, 6*levs+ 1)            surface pressure
c            (:, 6*levs+ 2)            surface geopotential
c            (:, 6*levs+ 3)            d(lnps)/dx
c            (:, 6*levs+ 4)            d(lnps)/dy
c            (:, 6*levs+ 5)            d(orog)/dx
c            (:, 6*levs+ 6)            d(orog)/dy
c
c            (:, 6*levs+ 7: 7*levs+ 6) geopotential energy
c            (:, 7*levs+ 7: 8*levs+ 6) x-gradient of tv
c            (:, 8*levs+ 7: 9*levs+ 6) y-gradient of tv
c            (:, 9*levs+ 7:10*levs+ 6) x-gradient of q
c            (:,10*levs+ 7:11*levs+ 6) y-gradient of q
c            (:,11*levs+ 7:12*levs+ 6) x-derivative of gz
c            (:,12*levs+ 7:13*levs+ 6) y-derivative of gz
c            (:,13*levs+ 7:13*levs+ntrace*levs+ 6) tracers
c            (:,13*levs+ 7:13*levs+ntrace*levs+ 6+ncloud*levs) clouds
c
c                                                                               
c   input files:                                                                
c     nss      - sigma spectral file                                            
c                                                                               
c subprograms called:                                                           
c   elat         compute latitudes                                              
c   fftfax       compute utility fields for fft                                 
c   gspc         compute utility fields for spectral transform                  
c   gradq        compute gradient in spectral space                             
c   dz2uv        compute vector components in spectral space                    
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
c
      dimension zz((jcap+1)*(jcap+2))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute utility fields                                                       
c                                                                               
      if(idrt.eq.0) then                                                        
        call elat(johf,slat,clat,wlat)                                          
      elseif(idrt.eq.4) then                                                    
        call glat(johf,slat,clat,wlat)                                          
      else                                                                      
        write(6,*) 'idrt should be 0 or 4 but is ',idrt                         
        call abort                                                              
      endif                                                                     
c                                                                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
            call fftfax (io2/2,ifax,trig)                                       
#endif
#ifdef DEFAULT
      call    fax(ifax,io2/2,3)                                                 
      call fftrig(trig,io2/2,3)                                                 
#endif
c                                                                               
      call gspc(jcap,eps,epstop,enn1,elonn1,eon,eontop)                         
c                                                                               
c  read sigma spectral data                                                     
c                                                                               
      nr=(jcap+1)*(jcap+2)                                                      
      read(nss) (ss(i,6*levs+2),i=1,nr)                                         
      read(nss) (ss(i,6*levs+1),i=1,nr)                                         
      do k=1,levs                                                               
        read(nss) (ss(i,2*levs+k),i=1,nr)                                       
      enddo                                                                     
      do k=1,levs                                                               
        read(nss) (ss(i,levs+k),i=1,nr)                                         
        read(nss) (ss(i,k),i=1,nr)                                              
      enddo                                                                     
      do k=1,levs                                                               
        read(nss) (ss(i,3*levs+k),i=1,nr)                                       
      enddo                                                                     
      if (ntrace+ncloud.ge.2.and.mflds.gt.18*levs+12) then
         do k=1,(ntrace+ncloud-1)*levs
            read(nss) (ss(i,18*levs+12+k),i=1,nr)
            write(6,*) 'rdss ss=',ss(1,18*levs+12+k)
         enddo
      endif
c
      do k=1,mflds                                                              
        do l=0,jcap                                                             
          sstop(2*l+1,k)=0.                                                     
          sstop(2*l+2,k)=0.                                                     
        enddo                                                                   
      enddo                                                                     
c
c  compute spectral coeff of geopotential height
c
      call hydrosp(nc,levs,si,sl,ss(1,6*levs+2),ss(1,2*levs+1),
     1             ss(1,6*levs+7))
c                                                                               
c  compute gradients and winds                                                  
c                                                                               
      call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,6*levs+2),                    
     &           ss(1,6*levs+5),ss(1,6*levs+6),sstop(1,6*levs+6))               
      call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,6*levs+1),                    
     &           ss(1,6*levs+3),ss(1,6*levs+4),sstop(1,6*levs+4))               
c     
      do k=1,levs                                                               
        call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,2*levs+k),
     &       ss(1, 7*levs+6+k),ss(1, 8*levs+6+k),sstop(1, 8*levs+6+k))
        call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,3*levs+k),
     &       ss(1, 9*levs+6+k),ss(1,10*levs+6+k),sstop(1,10*levs+6+k))
        call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,6*levs+6+k),
     &       ss(1,11*levs+6+k),ss(1,12*levs+6+k),sstop(1,12*levs+6+k))
        call dz2uv(jcap,enn1,elonn1,eon,eontop,ss(1,levs+k),ss(1,k),            
     &       ss(1,4*levs+k),ss(1,5*levs+k),                           
     &       sstop(1,4*levs+k),sstop(1,5*levs+k))                     
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
