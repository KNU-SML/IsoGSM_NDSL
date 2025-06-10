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
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
#include <define.h>
      subroutine fftidim (a, b, lot, isign)                                     
#include <paramodel.h>
      save                                                                      
#include <paramter.h>
      parameter(idimt=idim*2,mwvp2=(mwave+1)*2)                                 
      parameter (lotmin=32,lotmax=64,ncpu=8)                                    
      dimension a( idim ,lot)                                                   
      dimension b( idim ,lot)                                                   
c                                                                               
      dimension ifax(20), trigs( idim ,2)                                       
c                                                                               
      dimension work( idim ,lotmax,4), al( idimt ,lotmax)                       
c                                                                               
       data itest/0/                                                            
        if (itest.eq.0)  then                                                   
            itest=1                                                             
            jump = idim *2                                                      
c
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
            call fftfax (idim,ifax,trigs)                                       
#endif
#ifdef ASLES
#undef DEFAULT
c            print *,'new fft initializaiton starts'
            call ldfrmfb(idim,0,0,0,0,0,ifax,trigs,0,ierr)
#endif
#ifdef DEFAULT
            call    fax (ifax,idim,3)                                           
            call fftrig (trigs,idim,3)                                          
#endif
#ifdef ASLES
            if (ierr .ge. 3000)  print 120
            if (ierr .ge. 3000)  stop
#else
            if (ifax(1) .eq. -99)  print 120                                    
            if (ifax(1) .eq. -99)  stop     
#endif                                    
120         format (' error in fft idim .   idim  not factorable. ')            
            print 140, idim                                                     
140         format (' fftfax called in fft idim .  lonf = ',i4)                 
      endif                                                                     
c                                                                               
c                                                                               
      if (isign .eq.  1)  then                                                  
c                                                                               
c  multiple fast fourier transform - synthesis.  isign=1                        
c  good for zonal wave number  80 .                                             
c                                                                               
c     dimension a( idim ,lot)                                                   
c                                                                               
c   input - lot sets of complex coefficients in                                 
c           a(1,j), a(2,j), ..., a( mwavep *2,j), j=1,...,lot.                  
c           a( mwavep *2+1,j), ..., a( idim ,j), j=1,...,lot are not set        
c           before call fft idim .                                              
c                                                                               
c  output - lot sets of grid values in                                          
c           a(1,j), a(2,j), ..., a( idim ,j), j=1,...,lot.                      
c                                                                               
      nlot=max0(lot/ncpu,lotmin)                                                
      nlot=min0(nlot    ,lotmax)                                                
c
      do 460 i=1,lot,nlot                                                       
      lots = min0(nlot, lot-i+1)                                                
      do 360 j=i,i+lots-1                                                       
      do 320 l=1, mwvp2                                                         
             al(l,j-i+1) = a(l,j)                                               
320   continue                                                                  
      do 340 l= mwvp2+1 , idimt                                                 
             al(l,j-i+1) = 0.0                                                  
340   continue                                                                  
360   continue                                                                  
c                                                                               
c     call fft for systhesis.                                              
c                                                                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(al,work,trigs,ifax,1,jump, idim ,lots,1)                   
#endif
#ifdef ASLES
#undef DEFAULT
c        print *,'new fft, fft starts'
        call ldfrmbf(idim,lots,al,1,jump,-1,ifax,trigs,work,ierr)
       do jasl=1,lots
        do iasl=1,idim+2
          al(1+(iasl-1)*1+(jasl-1)*jump)=
     &      al(1+(iasl-1)*1+(jasl-1)*jump)/dble(idim)
        enddo
       enddo
#endif
#ifdef DEFAULT
        call fft99m (al,work,trigs,ifax,1,jump, idim ,lots,1)
#endif

#ifdef ASLES
            if (ierr .ge. 3000)  print 150
            if (ierr .ge. 3000)  stop
150         format (' error in asles FFT ')
#endif

c                                                                               
      do 440 j=i,i+lots-1                                                       
      do 420 l=1, idim                                                          
             a(l,j) = al(l,j-i+1)                                               
420   continue                                                                  
440   continue                                                                  
460   continue                                                                  
c                                                                               
      do 461 j=1,lot                                                            
      do 461 i=1,idim                                                           
      b(i,j)=a(i,j)                                                             
  461 continue                                                                  
c                                                                               
      endif                                                                     
c                                                                               
c                                                                               
      if (isign .eq. -1)  then                                                  
c                                                                               
c  multiple fast fourier transform - analysis.  isign=-1                        
c  good for zonal wave number  80 .                                             
c                                                                               
c     dimension a( idim ,lot), b( idim ,lot)                                    
c                                                                               
c   input - lot sets of grid values in                                          
c           a(1,j), a(2,j), ..., a( idim ,j), j=1,...,lot.                      
c           a array is not changed by subroutine fft idim .                     
c                                                                               
c  output - lot sets of complex coefficients in                                 
c           b(1,j), b(2,j), ..., b( mwavep *2,j), j=1,...,lot.                  
c           b( mwavep *2+1,j), ..., b( idim ,j), j=1,...,lot are not set        
c                                                                               
      nlot=max0(lot/ncpu,lotmin)                                                
      nlot=min0(nlot    ,lotmax)                                                
c
      do 660 i=1,lot,nlot                                                       
      lots = min0(nlot, lot-i+1)                                                
      do 560 j=i,i+lots-1                                                       
      do 520 l=1, idim                                                          
             al(l,j-i+1) = a(l,j)                                               
520   continue                                                                  
      do 540 l= idim+1 , idimt                                                  
             al(l,j-i+1) = 0.0                                                  
540   continue                                                                  
560   continue                                                                  
c                                                                               
c     call fft for analysis.                                               
c                                                                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
      call rfftmlt(al,work,trigs,ifax,1,jump, idim ,lots,-1)                  
#endif
#ifdef ASLES
#undef DEFAULT
c        print *,'new fft another fft starts'
        call ldfrmbf(idim,lots,al,1,jump,1,ifax,trigs,work,ierr)
       do jasl=1,lots
        do iasl=1,idim+2
c        do iasl=1,idim
          al(1+(iasl-1)*1+(jasl-1)*jump)=
     &      al(1+(iasl-1)*1+(jasl-1)*jump)/dble(idim)
        enddo
       enddo
#endif
#ifdef DEFAULT
      call fft99m (al,work,trigs,ifax,1,jump, idim ,lots,-1)                    
#endif

#ifdef ASLES
            if (ierr .ge. 3000)  print 150
            if (ierr .ge. 3000)  stop
#endif
c                                                                               
      do 640 j=i,i+lots-1                                                       
      do 620 l=1, mwvp2                                                         
             b(l,j) = al(l,j-i+1)                                               
620   continue                                                                  
640   continue                                                                  
660   continue                                                                  
c                                                                               
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                           
