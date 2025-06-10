      subroutine spline(l,n,x,f,s,p,xp,fp,dp)                                   
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    spline      interpolate data using cubic splines              
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: interpolate cubic spline profiles to given points.                  
c   cubic splines are piecewise cubic polynomials fitting the data              
c   with continuous first and second derivatives at interior points             
c   and second derivatives set to zero at and beyond the end points.            
c   subprogram spcoef must be already called to compute 2nd derivatives.        
c   note the inner dimension of the data is the number of profiles.             
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call spline(l,n,x,f,s,p,xp,fp,dp)                                   
c                                                                               
c   input argument list:                                                        
c     l        - integer number of profiles                                     
c     n        - integer number of points in each profile                       
c     x        - real (n) monotonically increasing abscissa values              
c     f        - real (l,n) data values                                         
c     s        - real (l,n) 2nd derivative of f (from subprogram spcoef)        
c     p        - real (l) point number or 0 to calculate point number           
c     xp       - real (l) abscissa values to which to interpolate               
c                                                                               
c   output argument list:                                                       
c     p        - real (l) point number or                                       
c     fp       - real (l) interpolated data values                              
c     dp       - real (l) 1st derivative of f at xp                             
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension x(n),f(l,n),s(l,n),p(l),xp(l),fp(l),dp(l)                       
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  calculate point number if necessary                                          
      do i=1,l                                                                  
        if(p(i).le.0.) then                                                     
          k=1                                                                   
          dowhile(k.le.n.and.xp(i).gt.x(k))                                     
            k=k+1                                                               
          enddo                                                                 
          p(i)=k-0.5                                                            
        endif                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  extrapolate or interpolate cubic spline                                      
      do i=1,l                                                                  
        if(p(i).le.1.) then                                                     
          p(i)=0.                                                               
          dx=x(2)-x(1)                                                          
          dp(i)=(f(i,2)-f(i,1))/dx-dx*s(i,2)/6.                                 
          fp(i)=f(i,1)+(xp(i)-x(1))*dp(i)                                       
        elseif(p(i).gt.n) then                                                  
          p(i)=n+1                                                              
          dx=x(n)-x(n-1)                                                        
          dp(i)=(f(i,n)-f(i,n-1))/dx+dx*s(i,n-1)/6.                             
          fp(i)=f(i,n)+(xp(i)-x(n))*dp(i)                                       
        else                                                                    
          kd=p(i)                                                               
          ku=kd+1                                                               
          dx=x(ku)-x(kd)                                                        
          dd=xp(i)-x(kd)                                                        
          du=dx-dd                                                              
          p(i)=kd+dd/dx                                                         
          fu=f(i,ku)                                                            
          fd=f(i,kd)                                                            
          df=fu-fd                                                              
          su=s(i,ku)                                                            
          sd=s(i,kd)                                                            
          ds=su-sd                                                              
          dp(i)=(df+su*dd**2/2-sd*du**2/2-ds*dx**2/6)/dx                        
          fp(i)=fd+dd/dx*(df-du*(dd*ds+dx*(su+2*sd))/6)                         
        endif                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
