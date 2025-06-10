      subroutine uv2dz(m,enn1,elonn1,eon,eontop,u,v,utop,vtop,d,z)              
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    uv2dz       compute divergence and vorticity from winds        
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes the divergence and vorticity from wind components          
c           in spectral space. subprogram gspc should be called already.        
c           if l is the zonal wavenumber, n is the total wavenumber,            
c           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
c           then the divergence d is computed as                                
c             d(l,n)=i*l*a*u(l,n)                                               
c                    +eps(l,n+1)*n*a*v(l,n+1)-eps(l,n)*(n+1)*a*v(l,n-1)         
c           and the vorticity z is computed as                                  
c             z(l,n)=i*l*a*v(l,n)                                               
c                    -eps(l,n+1)*n*a*u(l,n+1)+eps(l,n)*(n+1)*a*u(l,n-1)         
c           where u is the zonal wind and v is the meridional wind.             
c           u and v are weighted by the secant of latitude.                     
c           extra terms are used over top of the spectral triangle.             
c           advantage is taken of the fact that eps(l,l)=0                      
c           in order to vectorize over the entire spectral triangle.            
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call uv2dz(m,enn1,elonn1,eon,eontop,u,v,utop,vtop,d,z)              
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
c     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
c     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
c     eontop   - real (m+1) epsilon/n*a over top                                
c     u        - real ((m+1)*(m+2)) zonal wind (over coslat)                    
c     v        - real ((m+1)*(m+2)) merid wind (over coslat)                    
c     utop     - real (2*(m+1)) zonal wind (over coslat) over top               
c     vtop     - real (2*(m+1)) merid wind (over coslat) over top               
c                                                                               
c   output argument list:                                                       
c     d        - real ((m+1)*(m+2)) divergence                                  
c     z        - real ((m+1)*(m+2)) vorticity                                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
      real eon((m+1)*(m+2)/2),eontop(m+1)                                       
      real u((m+1)*(m+2)),v((m+1)*(m+2)),utop(2*(m+1)),vtop(2*(m+1))            
      real d((m+1)*(m+2)),z((m+1)*(m+2))                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute terms from the spectral triangle                                     
      i=1                                                                       
      d(2*i-1)=0.                                                               
      d(2*i)=0.                                                                 
      z(2*i-1)=0.                                                               
      z(2*i)=0.                                                                 
      do i=2,(m+1)*(m+2)/2-1                                                    
        d(2*i-1)=-elonn1(i)*u(2*i)+eon(i+1)*v(2*i+1)-eon(i)*v(2*i-3)            
        d(2*i)=elonn1(i)*u(2*i-1)+eon(i+1)*v(2*i+2)-eon(i)*v(2*i-2)             
        z(2*i-1)=-elonn1(i)*v(2*i)-eon(i+1)*u(2*i+1)+eon(i)*u(2*i-3)            
        z(2*i)=elonn1(i)*v(2*i-1)-eon(i+1)*u(2*i+2)+eon(i)*u(2*i-2)             
      enddo                                                                     
      i=(m+1)*(m+2)/2                                                           
      d(2*i-1)=-elonn1(i)*u(2*i)-eon(i)*v(2*i-3)                                
      d(2*i)=elonn1(i)*u(2*i-1)-eon(i)*v(2*i-2)                                 
      z(2*i-1)=-elonn1(i)*v(2*i)+eon(i)*u(2*i-3)                                
      z(2*i)=elonn1(i)*v(2*i-1)+eon(i)*u(2*i-2)                                 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  compute terms from over top of the spectral triangle                         
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
      do l=0,m                                                                  
        i=l*(2*m+1-l)/2+m+1                                                     
        d(2*i-1)=d(2*i-1)+eontop(l+1)*vtop(2*l+1)                               
        d(2*i)=d(2*i)+eontop(l+1)*vtop(2*l+2)                                   
        z(2*i-1)=z(2*i-1)-eontop(l+1)*utop(2*l+1)                               
        z(2*i)=z(2*i)-eontop(l+1)*utop(2*l+2)                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  multiply by laplacian term                                                   
      do i=2,(m+1)*(m+2)/2                                                      
        d(2*i-1)=d(2*i-1)*enn1(i)                                               
        d(2*i)=d(2*i)*enn1(i)                                                   
        z(2*i-1)=z(2*i-1)*enn1(i)                                               
        z(2*i)=z(2*i)*enn1(i)                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
