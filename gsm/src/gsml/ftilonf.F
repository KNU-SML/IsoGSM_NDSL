#define IGRD1P igrd1p
#define IGRD igrd
#define IGRD1 igrd1
#define IGRD2 igrd2
#define LATG latg
#define LONF22 lonf22
#define LONF lonf
#define LONFP lonfp
      subroutine ftilonf (a, b, lot, lcapf, lonff, lat, isign)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
#ifdef REDUCE_GRID
#define LATS latg_/2
#define LATX lat
#else
#define LATS 1
#define LATX 1
#endif
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
#ifdef MP
#ifdef RSM
#define JUMPR jumpr
#define NVECT (igrd1p_)*2
#else
#define NVECT lonfp_*2
#endif
#else
#ifdef RSM
#define JUMPR jumpr
#define NVECT (igrd_+1)*2
#else
#define NVECT 64
#endif
#endif
      parameter (lotmin=NVECT,lotmax=NVECT,ncpu=1)                            
c                                                                               
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
      parameter (jumph=(lonf_+2)/2,jump=jumph*2)
      dimension table(42000)
#ifndef RSM
      common /comlonf/crscale(LATS),rcscale(LATS)
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
      parameter (jump=lonf_ + 3)
      common /comlonf/scale(LATS),trig(lonf_+15,LATS)                        
#endif
#ifdef FFTW
#undef DEFAULT
      parameter (jump=lonf_ + 3)
      common /comlonf/ scale(LATS)
      common /comlonfi/iplan_c_to_r(LATS),iplan_r_to_c(LATS)
      dimension ak(jump,lotmax)
#endif
#ifdef DEFAULT
      parameter (jump=lonf_ + 3)
      dimension work(lonf_,lotmax,2)
      common /comlonf/ifax(20,LATS),trigs(lonf_,2,LATS)                        
#endif
      dimension al(jump,lotmax)
      dimension a(lonf22_,lot)                                                  
      dimension b(lonf22_,lot)                                                  
c                                                                               
      jcapf=lcapf-1
      twoj1=2*(jcapf+1)
c
      if (isign .eq.  1)  then                                                  
c                                                                               
c  multiple fast fourier transform - synthesis.  isign=1                        
c  good for zonal wave number jcap_.                                            
c                                                                               
c     dimension a(lonf_,lot)                                                    
c                                                                               
c   input - lot sets of complex coefficients in                                 
c           a(1,j), a(2,j), ..., a(jcap1_*2,j), j=1,...,lot.                    
c           a(jcap1_*2+1,j), ..., a(lonf_,j), j=1,...,lot are not set           
c           before call fftlonf.                                                
c                                                                               
c  output - lot sets of grid values in                                          
c           a(1,j), a(2,j), ..., a(lonf_,j), j=1,...,lot.                       
c                                                                               
c
      nlot=max0(lot/ncpu,lotmin)                                                
      nlot=min0(nlot    ,lotmax)                                                
c
      do 460 i=1,lot,nlot                                                       
        lots = min0(nlot, lot-i+1)                                                
        do 360 j=i,i+lots-1                                                       
          jj=j-i+1
          k=(j+1)/2                                                                 
          loff=mod(j+1,2)*lonf_                                                     
#ifdef FFTW
          al(1,jj)=a(1+loff,k)
          do l=2,jcapf+1
            al(l        ,jj)=a(2*l-1+loff,k)
            al(lonff+2-l,jj)=a(2*l  +loff,k)
          enddo
          do l=jcapf+2,lonff-jcapf
            al(l,jj)=0.0
          enddo
          do l=lonff+1,jump
            al(l,jj)=0.0
          enddo
#else
          do 320 l=1,twoj1                                                         
            al(l,jj) = a(l+loff,k)                                          
320       continue                                                                  
          do 340 l=twoj1 + 1 , jump
            al(l,jj) = 0.0                                                  
340       continue                                                                  
#endif
360     continue                                                                  
c                                                                               
c     call  fft for systhesis.                                              
c                                                                               
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifdef RSM
      scale=1.0
      call dcrft(1,al,jumph,al,jump,lonf_,lots,
     1           -1,scale,table,22000,table(22001),20000)
      call dcrft(0,al,jumph,al,jump,lonf_,lots,
     1           -1,scale,table,22000,table(22001),20000)
#else
      call dcrft(1,al,jumph,al,jump,lonff,lots,
     1           -1,crscale(LATX),table,22000,table(22001),20000)
      call dcrft(0,al,jumph,al,jump,lonff,lots,
     1           -1,crscale(LATX),table,22000,table(22001),20000)
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
       call zdfftm1du(-1,lonff,lots,al,1,jump,trig(1,LATX))                    
#endif
#ifdef FFTW
#undef DEFAULT
        call rfftw_f77(iplan_c_to_r(LATX),lots,al,1,jump,ak,1,jump)
#endif
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(al,work,trigs(1,1,LATX),ifax(1,LATX),
     1                                  1,jump,lonff,lots,1)                    
#endif
#ifdef ASLES
#undef DEFAULT
c        print *,'new fft, fft starts'
        call ldfrmbf(lonff,lots,al,1,jump,-1,ifax(1,LATX),
     1                            trigs(1,1,LATX),work,ierr)
#endif
#ifdef DEFAULT
        call fft99m (al,work,trigs(1,1,LATX),ifax(1,LATX),
     1                                  1,jump,lonff,lots,1)                      
#endif
#ifdef ASLES
            if (ierr .ge. 3000)  print 150
            if (ierr .ge. 3000)  stop
150         format (' error in asles FFT ')
#endif
c
#ifdef FFTW
#define ALK ak
#else
#define ALK al
#endif
        do 440 j=i,i+lots-1                                                       
          k=(j+1)/2                                                                 
          loff=mod(j+1,2)*lonff                                                     
          do 420 l=1,lonff                                                          
            b(l+loff,k) = ALK(l,j-i+1)                                          
420       continue                                                                  
440     continue                                                                  
460   continue                                                                  
c
      endif                                                                     
c                                                                               
c                                                                               
      if (isign .eq. -1)  then                                                  
c                                                                               
c  multiple fast fourier transform - analysis.  isign=-1                        
c  good for zonal wave number jcap_.                                            
c                                                                               
c     dimension a(lonf_,lot), b(lonf_,lot)                                      
c                                                                               
c   input - lot sets of grid values in                                          
c           a(1,j), a(2,j), ..., a(lonf_,j), j=1,...,lot.                       
c           a array is not changed by sr fftlonf.                               
c                                                                               
c  output - lot sets of complex coefficients in                                 
c           b(1,j), b(2,j), ..., b(jcap1_*2,j), j=1,...,lot.                    
c           b(jcap1_*2+1,j), ..., b(lonf_,j), j=1,...,lot are not set.          
c                                                                               
      nlot=max0(lot/ncpu,lotmin)                                                
      nlot=min0(nlot    ,lotmax)                                                
      do 660 i=1,lot,nlot                                                       
        lots = min0(nlot, lot-i+1)                                                
        do 560 j=i,i+lots-1                                                       
          k=(j+1)/2                                                                 
          loff=mod(j+1,2)*lonff                                                     
          do 520 l=1,lonff                                                          
             al(l,j-i+1) = a(l+loff,k)                                          
520       continue                                                                  
560     continue                                                                  
c                                                                               
c     call fft for analysis.                                               
c
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifdef RSM
        scale=1./float(lonf_)
        call drcft(1,al,jump,al,jumph,lonf_,lots,
     1             1,scale,table,22000,table(22001),20000)
        call drcft(0,al,jump,al,jumph,lonf_,lots,
     1             1,scale,table,22000,table(22001),20000)
#else
        call drcft(1,al,jump,al,jumph,lonff,lots,
     1             1,rcscale(LATX),table,22000,table(22001),20000)
        call drcft(0,al,jump,al,jumph,lonff,lots,
     1             1,rcscale(LATX),table,22000,table(22001),20000)
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
        call dzfftm1du(1,lonff,lots,al,1,jump,trig(1,LATX))                    
        call dscalm1d(lonff,lots,scale(LATX),al,1,jump)
#endif
#ifdef FFTW
#undef DEFAULT
        call rfftw_f77(iplan_r_to_c(LATX),lots,al,1,jump,ak,1,jump)
#endif
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(al,work,trigs(1,1,LATX),ifax(1,LATX),
     1                                  1,jump,lonff,lots,-1)                   
#endif
#ifdef ASLES
#undef DEFAULT
c        print *,'new fft, fft starts'
        call ldfrmbf(lonff,lots,al,1,jump,1,ifax(1,LATX),
     1                              trigs(1,1,LATX),work,ierr)
       do jasl=1,lots
        do iasl=1,lonff+2
          al(1+(iasl-1)*1+(jasl-1)*jump)=
     &      al(1+(iasl-1)*1+(jasl-1)*jump)/dble(lonff)
        enddo
       enddo
#endif
#ifdef DEFAULT
        call fft99m (al,work,trigs(1,1,LATX),ifax(1,LATX),
     1                                  1,jump,lonff,lots,-1)                     
#endif
#ifdef ASLES
            if (ierr .ge. 3000)  print 160
            if (ierr .ge. 3000)  stop
160         format (' error in asles FFT ')
#endif
c                                                                               
        do 640 j=i,i+lots-1                                                       
          k=(j+1)/2                                                                 
          loff=mod(j+1,2)*lonf_                                                     
#ifdef FFTW
          jj=j-i+1
          b(1+loff,k)      =ak(1,jj)*scale
          b(2+loff,k)      =0.0
          do l=2,jcapf+1
            b(2*l-1+loff,k)=ak(l        ,jj)*scale
            b(2*l  +loff,k)=ak(lonff+2-l,jj)*scale
          enddo
#else
          do 620 l=1,twoj1                                                         
             b(l+loff,k) = al(l,j-i+1)                                          
620       continue                                                                  
#endif
640     continue                                                                  
660   continue                                                                  
c
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
