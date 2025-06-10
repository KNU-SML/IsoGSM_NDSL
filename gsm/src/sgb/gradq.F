      subroutine gradq(m,enn1,elonn1,eon,eontop,q,                              
     &                 qdx,qdy,qdytop)                                          
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    gradq       compute gradient in spectral space                 
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: computes the horizontal vector gradient of a scalar field           
c           in spectral space. subprogram gspc should be called already.        
c           if l is the zonal wavenumber, n is the total wavenumber,            
c           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
c           then the zonal gradient of q(l,n) is simply i*l/a*q(l,n)            
c           while the meridional gradient of q(l,n) is computed as              
c           eps(l,n+1)*(n+2)/a*q(l,n+1)-eps(l,n+1)*(n-1)/a*q(l,n-1).            
c           extra terms are computed over top of the spectral triangle.         
c           advantage is taken of the fact that eps(l,l)=0                      
c           in order to vectorize over the entire spectral triangle.            
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call gradq(m,enn1,elonn1,eon,eontop,q,                              
c    &                 qdx,qdy,qdytop)                                          
c                                                                               
c   input argument list:                                                        
c     m        - integer spectral truncation                                    
c     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
c     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
c     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
c     eontop   - real (m+1) epsilon/n*a over top                                
c     q        - real ((m+1)*(m+2)) scalar field                                
c                                                                               
c   output argument list:                                                       
c     qdx      - real ((m+1)*(m+2)) zonal gradient (times coslat)               
c     qdy      - real ((m+1)*(m+2)) merid gradient (times coslat)               
c     qdytop   - real (2*(m+1)) merid gradient (times coslat) over top          
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
      real eon((m+1)*(m+2)/2),eontop(m+1)                                       
      real q((m+1)*(m+2))                                                       
      real qdx((m+1)*(m+2)),qdy((m+1)*(m+2)),qdytop(2*(m+1))                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  take zonal and meridional gradients                                          
      i=1                                                                       
      qdx(2*i-1)=0.                                                             
      qdx(2*i)=0.                                                               
      qdy(2*i-1)=eon(i+1)*enn1(i+1)*q(2*i+1)                                    
      qdy(2*i)=eon(i+1)*enn1(i+1)*q(2*i+2)                                      
      do i=2,(m+1)*(m+2)/2-1                                                    
        qdx(2*i-1)=-elonn1(i)*enn1(i)*q(2*i)                                    
        qdx(2*i)=elonn1(i)*enn1(i)*q(2*i-1)                                     
        qdy(2*i-1)=eon(i+1)*enn1(i+1)*q(2*i+1)-eon(i)*enn1(i-1)*q(2*i-3)        
        qdy(2*i)=eon(i+1)*enn1(i+1)*q(2*i+2)-eon(i)*enn1(i-1)*q(2*i-2)          
      enddo                                                                     
      i=(m+1)*(m+2)/2                                                           
      qdx(2*i-1)=-elonn1(i)*enn1(i)*q(2*i)                                      
      qdx(2*i)=elonn1(i)*enn1(i)*q(2*i-1)                                       
      qdy(2*i-1)=-eon(i)*enn1(i-1)*q(2*i-3)                                     
      qdy(2*i)=-eon(i)*enn1(i-1)*q(2*i-2)                                       
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  take meridional gradient over top                                            
      do l=0,m                                                                  
        i=l*(2*m+1-l)/2+m+1                                                     
        qdytop(2*l+1)=-eontop(l+1)*enn1(i)*q(2*i-1)                             
        qdytop(2*l+2)=-eontop(l+1)*enn1(i)*q(2*i)                               
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
