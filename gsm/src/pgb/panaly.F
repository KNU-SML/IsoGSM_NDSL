#define KO ko
#include <define.h>
      subroutine panaly(m,im,ix,nc,nctop,km,wgt,clat,pln,plntop,mp,             
     &                  f,spc,spctop)                                           
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    panaly      analyze spectral from fourier                      
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: analyzes spectral coefficients from fourier coefficients            
c           for a latitude pair (northern and southern hemispheres).            
c           vector components are multiplied by cosine of latitude.             
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c   94-08-01  mark iredell   moved zonal wavenumber loop inside                 
c                                                                               
c usage:    call panaly(m,im,ix,nc,nctop,km,wgt,clat,pln,plntop,mp,             
c    &                  f,spc,spctop)                                           
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
c     wgt      - real gaussian weight                                           
c     clat     - real cosine of latitude                                        
c     pln      - real ((m+1)*(m+2)/2) legendre polynomials                      
c     plntop   - real (m+1) legendre polynomial over top                        
c     mp       - integer (km) identifiers (0 for scalar, 1 for vector)          
c     f        - real (ix,2,km) fourier coefficients combined                   
c     spc      - real (nc,km) spectral coefficients                             
c     spctop   - real (nctop,km) spectral coefficients over top                 
c                                                                               
c   output argument list:                                                       
c     spc      - real (nc,km) spectral coefficients                             
c     spctop   - real (nctop,km) spectral coefficients over top                 
c                                                                               
c subprograms called:                                                           
c   sgerx1       cray library matrix rank 1 update                              
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      integer mp(km)                                                            
      real pln((m+1)*(m+2)/2),plntop(m+1)                                       
      real f(ix,2,km)                                                           
      real spc(nc,km),spctop(nctop,km)                                          
c                                                                               
#define DEFAULT
#ifdef DYNAMIC_ALLOC
#undef DEFAULT
      real fw(2,2,km)                                                           
#endif
#ifdef DEFAULT
      real fw(2,2,ko_*5)                                                        
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  for each zonal wavenumber, analyze terms over total wavenumber.              
c  analyze even and odd polynomials separately.                                 
c  commented code replaced by library calls.                                    
      lx=min(m,im/2)                                                            
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 shared(lx,m,km,mp,wgt,clat,f,spctop,plntop,spc,pln,nc)                   
cmic$1 private(l,nt,k,fw,is,ip,n,k)                                             
#endif
      do l=0,lx                                                                 
        nt=mod(m+1-l,2)+1                                                       
        do k=1,km                                                               
          if(mp(k).eq.0) then                                                   
            fw(1,1,k)=wgt*(f(2*l+1,1,k)+f(2*l+1,2,k))                           
            fw(2,1,k)=wgt*(f(2*l+2,1,k)+f(2*l+2,2,k))                           
            fw(1,2,k)=wgt*(f(2*l+1,1,k)-f(2*l+1,2,k))                           
            fw(2,2,k)=wgt*(f(2*l+2,1,k)-f(2*l+2,2,k))                           
          else                                                                  
            fw(1,1,k)=wgt*clat*(f(2*l+1,1,k)+f(2*l+1,2,k))                      
            fw(2,1,k)=wgt*clat*(f(2*l+2,1,k)+f(2*l+2,2,k))                      
            fw(1,2,k)=wgt*clat*(f(2*l+1,1,k)-f(2*l+1,2,k))                      
            fw(2,2,k)=wgt*clat*(f(2*l+2,1,k)-f(2*l+2,2,k))                      
            spctop(2*l+1,k)=spctop(2*l+1,k)+plntop(l+1)*fw(1,nt,k)              
            spctop(2*l+2,k)=spctop(2*l+2,k)+plntop(l+1)*fw(2,nt,k)              
          endif                                                                 
        enddo                                                                   
        is=l*(2*m+1-l)                                                          
        ip=is/2+1                                                               
#define DEFAULT
#ifdef SGERX1
#undef DEFAULT
        call sgerx1((m+2-l)/2,km,1.,pln(ip+l),2,fw(1,1,1),4,                    
     &              spc(is+2*l+1,1),4,nc)                                       
        call sgerx1((m+2-l)/2,km,1.,pln(ip+l),2,fw(2,1,1),4,                    
     &              spc(is+2*l+2,1),4,nc)                                       
        call sgerx1((m+1-l)/2,km,1.,pln(ip+l+1),2,fw(1,2,1),4,                  
     &              spc(is+2*l+3,1),4,nc)                                       
        call sgerx1((m+1-l)/2,km,1.,pln(ip+l+1),2,fw(2,2,1),4,                  
     &              spc(is+2*l+4,1),4,nc)                                       
#endif
#ifdef DEFAULT
        do n=l,m,2                                                              
          do k=1,km                                                             
            spc(is+2*n+1,k)=spc(is+2*n+1,k)+pln(ip+n)*fw(1,1,k)                 
            spc(is+2*n+2,k)=spc(is+2*n+2,k)+pln(ip+n)*fw(2,1,k)                 
          enddo                                                                 
        enddo                                                                   
        do n=l+1,m,2                                                            
          do k=1,km                                                             
            spc(is+2*n+1,k)=spc(is+2*n+1,k)+pln(ip+n)*fw(1,2,k)                 
            spc(is+2*n+2,k)=spc(is+2*n+2,k)+pln(ip+n)*fw(2,2,k)                 
          enddo                                                                 
        enddo                                                                   
#endif
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
