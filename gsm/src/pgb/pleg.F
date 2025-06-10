      subroutine pleg(m,slat,clat,eps,epstop,pln,plntop)                        
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    pleg        compute legendre polynomials                       
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: evaluates the orthonormal associated legendre polynomials           
c           in the spectral triangle at a given latitude.                       
c           subprogram gspc should be called already.                           
c           if l is the zonal wavenumber, n is the total wavenumber,            
c           and eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) then                      
c           the following bootstrapping formulas are used:                      
c           pln(0,0)=sqrt(0.5)                                                  
c           pln(l,l)=pln(l-1,l-1)*clat*sqrt(float(2*l+1)/float(2*l))            
c           pln(l,n)=(slat*pln(l,n-1)-eps(l,n-1)*pln(l,n-2))/eps(l,n)           
c           synthesis at the pole needs only two zonal wavenumbers.             
c           scalar fields are synthesized with zonal wavenumber 0 while         
c           vector fields are synthesized with zonal wavenumber 1.              
c           (thus polar vector fields are implicitly divided by clat.)          
c           the following bootstrapping formulas are used at the pole:          
c           pln(0,0)=sqrt(0.5)                                                  
c           pln(1,1)=sqrt(0.75)                                                 
c           pln(l,n)=(pln(l,n-1)-eps(l,n-1)*pln(l,n-2))/eps(l,n)                
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call pleg(m,slat,clat,eps,epstop,pln,plntop)                        
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c     slat     - real sine of latitude                                          
c     clat     - real cosine of latitude                                        
c     eps      - real ((m+1)*(m+2)/2) sqrt((n**2-l**2)/(4*n**2-1))              
c     epstop   - real (m+1) sqrt((n**2-l**2)/(4*n**2-1)) over top               
c                                                                               
c   output argument list:                                                       
c     pln      - real ((m+1)*(m+2)/2) legendre polynomial                       
c     plntop   - real (m+1) legendre polynomial over top                        
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
      real eps((m+1)*(m+2)/2),epstop(m+1)                                       
      real pln((m+1)*(m+2)/2),plntop(m+1)                                       
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  iteratively compute pln within spectral triangle at pole                     
      if(clat.eq.0.) then                                                       
        pln(1)=sqrt(0.5)                                                        
        pln(m+2)=sqrt(0.75)                                                     
        pln(2)=pln(1)/eps(2)                                                    
        pln(m+3)=pln(m+2)/eps(m+3)                                              
        pln(3)=(pln(2)-eps(2)*pln(1))/eps(3)                                    
        do n=3,m                                                                
          i=n+1                                                                 
          pln(i)=(pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                            
          i=n+m+1                                                               
          pln(i)=(pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                            
        enddo                                                                   
        do i=2*m+2,(m+1)*(m+2)/2                                                
          pln(i)=0.                                                             
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute polynomials over top of spectral triangle                            
        i=m+2                                                                   
        plntop(1)=(pln(i-1)-eps(i-1)*pln(i-2))/epstop(1)                        
        i=2*m+2                                                                 
        plntop(2)=(pln(i-1)-eps(i-1)*pln(i-2))/epstop(2)                        
        do l=2,m                                                                
          plntop(l+1)=0.                                                        
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  iteratively compute pln(l,l) (bottom hypotenuse of triangle)                 
      else                                                                      
        nml=0                                                                   
        i=1                                                                     
        pln(i)=sqrt(0.5)                                                        
        do l=1,m-nml                                                            
          plni=pln(i)                                                           
          i=l*(2*m+3-l)/2+(nml+1)                                               
          pln(i)=plni*clat*sqrt(float(2*l+1)/float(2*l))                        
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute pln(l,l+1) (diagonal next to bottom hypotenuse of triangle)          
        nml=1                                                                   
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
        do l=0,m-nml                                                            
          i=l*(2*m+3-l)/2+(nml+1)                                               
          pln(i)=slat*pln(i-1)/eps(i)                                           
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute remaining pln in spectral triangle                                   
        do nml=2,m                                                              
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
          do l=0,m-nml                                                          
            i=l*(2*m+3-l)/2+(nml+1)                                             
            pln(i)=(slat*pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                     
          enddo                                                                 
        enddo                                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute polynomials over top of spectral triangle                            
        do l=0,m                                                                
          nml=m+1-l                                                             
          i=l*(2*m+3-l)/2+(nml+1)                                               
          plntop(l+1)=(slat*pln(i-1)-eps(i-1)*pln(i-2))/epstop(l+1)             
        enddo                                                                   
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
