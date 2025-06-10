      subroutine dz2uv(m,enn1,elonn1,eon,eontop,d,z,                            
     &                 u,v,utop,vtop)                                           
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    dz2uv       compute winds from divergence and vorticity        
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes the wind components from divergence and vorticity          
c           in spectral space. subprogram gspc should be called already.        
c           if l is the zonal wavenumber, n is the total wavenumber,            
c           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
c           then the zonal wind component u is computed as                      
c             u(l,n)=-i*l/(n*(n+1))*a*d(l,n)                                    
c                    +eps(l,n+1)/(n+1)*a*z(l,n+1)-eps(l,n)/n*a*z(l,n-1)         
c           and the meridional wind component v is computed as                  
c             v(l,n)=-i*l/(n*(n+1))*a*z(l,n)                                    
c                    -eps(l,n+1)/(n+1)*a*d(l,n+1)+eps(l,n)/n*a*d(l,n-1)         
c           where d is divergence and z is vorticity.                           
c           extra terms are computed over top of the spectral triangle.         
c           advantage is taken of the fact that eps(l,l)=0                      
c           in order to vectorize over the entire spectral triangle.            
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call dz2uv(m,enn1,elonn1,eon,eontop,d,z,                            
c    &                 u,v,utop,vtop)                                           
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
c     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
c     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
c     eontop   - real (m+1) epsilon/n*a over top                                
c     d        - real ((m+1)*(m+2)) divergence                                  
c     z        - real ((m+1)*(m+2)) vorticity                                   
c                                                                               
c   output argument list:                                                       
c     u        - real ((m+1)*(m+2)) zonal wind (times coslat)                   
c     v        - real ((m+1)*(m+2)) merid wind (times coslat)                   
c     utop     - real (2*(m+1)) zonal wind (times coslat) over top              
c     vtop     - real (2*(m+1)) merid wind (times coslat) over top              
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
      real eon((m+1)*(m+2)/2),eontop(m+1)                                       
      real d((m+1)*(m+2)),z((m+1)*(m+2))                                        
      real u((m+1)*(m+2)),v((m+1)*(m+2)),utop(2*(m+1)),vtop(2*(m+1))            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute winds in the spectral triangle                                       
      i=1                                                                       
      u(2*i-1)=eon(i+1)*z(2*i+1)                                                
      u(2*i)=eon(i+1)*z(2*i+2)                                                  
      v(2*i-1)=-eon(i+1)*d(2*i+1)                                               
      v(2*i)=-eon(i+1)*d(2*i+2)                                                 
      do i=2,(m+1)*(m+2)/2-1                                                    
        u(2*i-1)=elonn1(i)*d(2*i)+eon(i+1)*z(2*i+1)-eon(i)*z(2*i-3)             
        u(2*i)=-elonn1(i)*d(2*i-1)+eon(i+1)*z(2*i+2)-eon(i)*z(2*i-2)            
        v(2*i-1)=elonn1(i)*z(2*i)-eon(i+1)*d(2*i+1)+eon(i)*d(2*i-3)             
        v(2*i)=-elonn1(i)*z(2*i-1)-eon(i+1)*d(2*i+2)+eon(i)*d(2*i-2)            
      enddo                                                                     
      i=(m+1)*(m+2)/2                                                           
      u(2*i-1)=elonn1(i)*d(2*i)-eon(i)*z(2*i-3)                                 
      u(2*i)=-elonn1(i)*d(2*i-1)-eon(i)*z(2*i-2)                                
      v(2*i-1)=elonn1(i)*z(2*i)+eon(i)*d(2*i-3)                                 
      v(2*i)=-elonn1(i)*z(2*i-1)+eon(i)*d(2*i-2)                                
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute winds over top of the spectral triangle                              
      do l=0,m                                                                  
        i=l*(2*m+1-l)/2+m+1                                                     
        utop(2*l+1)=-eontop(l+1)*z(2*i-1)                                       
        utop(2*l+2)=-eontop(l+1)*z(2*i)                                         
        vtop(2*l+1)=eontop(l+1)*d(2*i-1)                                        
        vtop(2*l+2)=eontop(l+1)*d(2*i)                                          
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
