#define IO io
#include <define.h>
      subroutine psynthv(m,im,ix,nc,nctop,km,clat,pln,plntop,mp,                
     &                  spc,spctop,f)                                           
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    psynth      synthesize fourier from spectral                   
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: synthesizes fourier coefficients from spectral coefficients         
c           for a latitude pair (northern and southern hemispheres).            
c           vector components are divided by cosine of latitude.                
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call psynth(m,im,ix,nc,nctop,km,clat,pln,plntop,mp,                 
c    &                  spc,spctop,f)                                           
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c     im       - integer even number of fourier coefficients                    
c     ix       - integer dimension of fourier coefficients (ix>=im+2)           
c     nc       - integer dimension of spectral coefficients                     
c                (nc>=(m+1)*(m+2))                                              
c     nctop    - integer dimension of spectral coefficients over top            
c                (nctop>=2*(m+1))                                               
c     km       - integer number of fields                                       
c     clat     - real cosine of latitude                                        
c     pln      - real ((m+1)*(m+2)/2) legendre polynomial                       
c     plntop   - real (m+1) legendre polynomial over top                        
c     spc      - real (nc,km) spectral coefficients                             
c     spctop   - real (nctop,km) spectral coefficients over top                 
c     mp       - integer (km) identifiers (0 for scalar, 1 for vector)          
c                                                                               
c   output argument list:                                                       
c     f        - real (ix,2,km) fourier coefficients for latitude pair          
c                                                                               
c subprograms called:                                                           
c   sgemvx1      cray library matrix times vector                               
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
      real pln((m+1)*(m+2)/2),plntop(m+1)                                       
      integer mp(km)                                                            
      real spc(nc,km),spctop(nctop,km)                                          
      real f(ix,2,km)                                                           
      mstop = min(io_/2, m)                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  synthesis over pole.                                                         
c  zero out fourier waves.                                                      
      if(clat.eq.0) then                                                        
        do k=1,km                                                               
          do l=0,im/2                                                           
            f(2*l+1,1,k)=0.                                                     
            f(2*l+2,1,k)=0.                                                     
            f(2*l+1,2,k)=0.                                                     
            f(2*l+2,2,k)=0.                                                     
          enddo                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  initialize fourier coefficients with terms over top of the spectrum.         
c  initialize even and odd polynomials separately.                              
        ltope=mod(m+1,2)                                                        
        do k=1,km                                                               
          l=mp(k)                                                               
          if(l.eq.ltope) then                                                   
            f(2*l+1,1,k)=plntop(l+1)*spctop(2*l+1,k)                            
            f(2*l+2,1,k)=plntop(l+1)*spctop(2*l+2,k)                            
            f(2*l+1,2,k)=0.                                                     
            f(2*l+2,2,k)=0.                                                     
          else                                                                  
            f(2*l+1,1,k)=0.                                                     
            f(2*l+2,1,k)=0.                                                     
            f(2*l+1,2,k)=plntop(l+1)*spctop(2*l+1,k)                            
            f(2*l+2,2,k)=plntop(l+1)*spctop(2*l+2,k)                            
          endif                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  for each zonal wavenumber, synthesize terms over total wavenumber.           
c  synthesize even and odd polynomials separately.                              
        do k=1,km                                                               
          l=mp(k)                                                               
          is=l*(2*m+1-l)                                                        
          ip=is/2+1                                                             
          do n=l,mstop,2                                                        
            f(2*l+1,1,k)=f(2*l+1,1,k)+pln(ip+n)*spc(is+2*n+1,k)                 
            f(2*l+2,1,k)=f(2*l+2,1,k)+pln(ip+n)*spc(is+2*n+2,k)                 
          enddo                                                                 
          do n=l+1,mstop,2                                                      
            f(2*l+1,2,k)=f(2*l+1,2,k)+pln(ip+n)*spc(is+2*n+1,k)                 
            f(2*l+2,2,k)=f(2*l+2,2,k)+pln(ip+n)*spc(is+2*n+2,k)                 
          enddo                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  separate fourier coefficients from each hemisphere.                          
c  odd polynomials contribute negatively to the southern hemisphere.            
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
        do k=1,km                                                               
          l=mp(k)                                                               
          f1r=f(2*l+1,1,k)                                                      
          f1i=f(2*l+2,1,k)                                                      
          f(2*l+1,1,k)=f1r+f(2*l+1,2,k)                                         
          f(2*l+2,1,k)=f1i+f(2*l+2,2,k)                                         
          f(2*l+1,2,k)=f1r-f(2*l+1,2,k)                                         
          f(2*l+2,2,k)=f1i-f(2*l+2,2,k)                                         
        enddo                                                                   
c                                                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  synthesis over finite latitude.                                              
c  initialize fourier coefficients with terms over top of the spectrum.         
c  initialize even and odd polynomials separately.                              
c                                                                               
      else                                                                      
        lx=min(m,im/2)                                                          
        ltope=mod(m+1,2)                                                        
        ltopo=1-ltope                                                           
        do k=1,km                                                               
          if(mp(k).eq.0) then                                                   
            do l=ltope,lx,2                                                     
              f(2*l+1,1,k)=0.                                                   
              f(2*l+2,1,k)=0.                                                   
              f(2*l+1,2,k)=0.                                                   
              f(2*l+2,2,k)=0.                                                   
            enddo                                                               
            do l=ltopo,lx,2                                                     
              f(2*l+1,1,k)=0.                                                   
              f(2*l+2,1,k)=0.                                                   
              f(2*l+1,2,k)=0.                                                   
              f(2*l+2,2,k)=0.                                                   
            enddo                                                               
          else                                                                  
            do l=ltope,lx,2                                                     
              f(2*l+1,1,k)=plntop(l+1)*spctop(2*l+1,k)                          
              f(2*l+2,1,k)=plntop(l+1)*spctop(2*l+2,k)                          
              f(2*l+1,2,k)=0.                                                   
              f(2*l+2,2,k)=0.                                                   
            enddo                                                               
            do l=ltopo,lx,2                                                     
              f(2*l+1,1,k)=0.                                                   
              f(2*l+2,1,k)=0.                                                   
              f(2*l+1,2,k)=plntop(l+1)*spctop(2*l+1,k)                          
              f(2*l+2,2,k)=plntop(l+1)*spctop(2*l+2,k)                          
            enddo                                                               
          endif                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  for each zonal wavenumber, synthesize terms over total wavenumber.           
c  synthesize even and odd polynomials separately.                              
c  commented code replaced by library calls.                                    
        do l=0,lx                                                               
          is=l*(2*m+1-l)                                                        
          ip=is/2+1                                                             
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
          call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+1,1),nc,4,                    
     &                 pln(ip+l),2,1.,f(2*l+1,1,1),ix*2)                        
          call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+2,1),nc,4,                    
     &                 pln(ip+l),2,1.,f(2*l+2,1,1),ix*2)                        
          call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+3,1),nc,4,                    
     &                 pln(ip+l+1),2,1.,f(2*l+1,2,1),ix*2)                      
          call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+4,1),nc,4,                    
     &                 pln(ip+l+1),2,1.,f(2*l+2,2,1),ix*2)                      
#endif
#ifdef DEFAULT
          do n=l,mstop,2                                                        
            do k=1,km                                                           
              f(2*l+1,1,k)=f(2*l+1,1,k)+pln(ip+n)*spc(is+2*n+1,k)               
              f(2*l+2,1,k)=f(2*l+2,1,k)+pln(ip+n)*spc(is+2*n+2,k)               
            enddo                                                               
          enddo                                                                 
            do n=l+1,mstop,2
            do k=1,km                                                           
              f(2*l+1,2,k)=f(2*l+1,2,k)+pln(ip+n)*spc(is+2*n+1,k)               
              f(2*l+2,2,k)=f(2*l+2,2,k)+pln(ip+n)*spc(is+2*n+2,k)               
            enddo                                                               
          enddo                                                                 
#endif
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  separate fourier coefficients from each hemisphere.                          
c  odd polynomials contribute negatively to the southern hemisphere.            
c  divide vector components by cosine latitude.                                 
        do k=1,km                                                               
          do l=0,lx                                                             
            f1r=f(2*l+1,1,k)                                                    
            f1i=f(2*l+2,1,k)                                                    
            f(2*l+1,1,k)=f1r+f(2*l+1,2,k)                                       
            f(2*l+2,1,k)=f1i+f(2*l+2,2,k)                                       
            f(2*l+1,2,k)=f1r-f(2*l+1,2,k)                                       
            f(2*l+2,2,k)=f1i-f(2*l+2,2,k)                                       
          enddo                                                                 
          if(mp(k).eq.1) then                                                   
            do l=0,lx                                                           
              f(2*l+1,1,k)=f(2*l+1,1,k)/clat                                    
              f(2*l+2,1,k)=f(2*l+2,1,k)/clat                                    
              f(2*l+1,2,k)=f(2*l+1,2,k)/clat                                    
              f(2*l+2,2,k)=f(2*l+2,2,k)/clat                                    
            enddo                                                               
          endif                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  zero out fourier waves outside of spectrum                                   
        do l=lx+1,im/2                                                          
          do k=1,km                                                             
            f(2*l+1,1,k)=0.                                                     
            f(2*l+2,1,k)=0.                                                     
            f(2*l+1,2,k)=0.                                                     
            f(2*l+2,2,k)=0.                                                     
          enddo                                                                 
        enddo                                                                   
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
