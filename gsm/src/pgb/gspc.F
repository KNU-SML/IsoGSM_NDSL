      subroutine gspc(m,eps,epstop,enn1,elonn1,eon,eontop)                      
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    gspc        compute utility spectral fields                    
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes constant fields indexed in the spectral triangle           
c           in "ibm order" (zonal wavenumber is the slower index).              
c           if l is the zonal wavenumber and n is the total wavenumber          
c           and a is the earth radius, then the fields returned are:            
c           (1) normalizing factor epsilon=sqrt((n**2-l**2)/(4*n**2-1))         
c           (2) laplacian factor n*(n+1)/a**2                                   
c           (3) zonal derivative/laplacian factor l/(n*(n+1))*a                 
c           (4) meridional derivative/laplacian factor epsilon/n*a              
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call gspc(m,eps,epstop,enn1,elonn1,eon,eontop)                      
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c                                                                               
c   output argument list:                                                       
c     eps      - real ((m+1)*(m+2)/2) sqrt((n**2-l**2)/(4*n**2-1))              
c     epstop   - real (m+1) sqrt((n**2-l**2)/(4*n**2-1)) over top               
c     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
c     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
c     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
c     eontop   - real (m+1) epsilon/n*a over top                                
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      real eps((m+1)*(m+2)/2),epstop(m+1)                                       
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
      real eon((m+1)*(m+2)/2),eontop(m+1)                                       
      parameter(rerth=6.3712e6,ra2=1./rerth**2)                                 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do l=0,m                                                                  
        ill=l*(2*m+3-l)/2+1                                                     
        eps(ill)=0.                                                             
        enn1(ill)=ra2*l*(l+1)                                                   
        elonn1(ill)=rerth/(l+1)                                                 
        eon(ill)=0.                                                             
      enddo                                                                     
      do l=0,m                                                                  
        is=l*(2*m+1-l)                                                          
        ip=is/2+1                                                               
        do n=l+1,m                                                              
          eps(ip+n)=sqrt(float(n**2-l**2)/float(4*n**2-1))                      
          enn1(ip+n)=ra2*n*(n+1)                                                
          elonn1(ip+n)=rerth*l/(n*(n+1))                                        
          eon(ip+n)=rerth/n*eps(ip+n)                                           
        enddo                                                                   
      enddo                                                                     
      do l=0,m                                                                  
        epstop(l+1)=sqrt(float((m+1)**2-l**2)/float(4*(m+1)**2-1))              
        eontop(l+1)=rerth/(m+1)*epstop(l+1)                                     
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
