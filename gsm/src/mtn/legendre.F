      subroutine legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)                     
c subprogram documentation block                                            
c                                                                               
c subprogram:  legendre   compute legendre polynomials.                         
c   prgmmr: iredell       org: w/nmc23       date: 92-04-16                     
c                                                                               
c abstract: evaluates the normalized associated legendre polynomials            
c   at a given cosine of colatitude.                                            
c                                                                               
c program history log:                                                          
c   92-04-16  iredell                                                           
c                                                                               
c usage:    call legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)                     
c                                                                               
c   input argument list:                                                        
c     ipd      - ipd = 0 to evaluate polynomial itself,                         
c                ipd = 1 to evaluate polynomial derivative.                     
c     cosclt   - cosine of colatitude at which to evalutate.                    
c     maxwv    - spectral truncation                                            
c     iromb    - iromb = 0 for triangular truncation                            
c                iromb = 1 for rhomboidal truncation                            
c                                                                               
c   work argument list:                                                         
c     ex       - real (maxwv+2) work area                                       
c     px       - real (maxwv+3) work area                                       
c                                                                               
c   output argument list:                                                       
c     pnm      - real (kmax) polynomial values,                                 
c                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
c                                                                               
c   remarks: fortran 9x extensions are used.                                    
c                                                                               
c attributes:                                                                   
c   cray ymp.                                                                   
c                                                                               
c
      save                                                                      
      dimension pnm((maxwv+1)*((iromb+1)*maxwv+2)/2)                            
      dimension ex(0:maxwv+1),px(-1:maxwv+1)                                    
      sinclt=sqrt(1.-cosclt**2)                                                 
      k=0                                                                       
      ex(0)=0.                                                                  
      px(-1)=0.                                                                 
      px(0)=sqrt(0.5)                                                           
      do m=0,maxwv                                                              
        if(m.gt.0) px(0)=px(0)*sinclt/(ex(1)*sqrt(float(2*m)))                  
        do n=m+1,maxwv+iromb*m+1                                                
          ex(n-m)=sqrt(float(n**2-m**2)/float(4*n**2-1))                        
        enddo                                                                   
        do n=m+1,maxwv+iromb*m+1                                                
          px(n-m)=(cosclt*px(n-m-1)-ex(n-m-1)*px(n-m-2))/ex(n-m)                
        enddo                                                                   
        if(ipd.eq.0) then                                                       
          do n=m,maxwv+iromb*m                                                  
            k=k+1                                                               
            pnm(k)=px(n-m)                                                      
          enddo                                                                 
        else                                                                    
          do n=m,maxwv+iromb*m                                                  
            k=k+1                                                               
            pnm(k)=n*ex(n-m+1)*px(n-m+1)-(n+1)*ex(n-m)*px(n-m-1)                
          enddo                                                                 
        endif                                                                   
      enddo                                                                     
      return                                                                    
      end                                                                       
