      subroutine spfmax(l,n,x,f,s,p,xp,fp)                                      
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    spfmax      find maximum value using cubic splines            
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: compute the maximum data value of cubic spline profiles.            
c   cubic splines are piecewise cubic polynomials fitting the data              
c   with continuous first and second derivatives at interior points             
c   and second derivatives set to zero at and beyond the end points.            
c   subprogram spcoef must be already called to compute 2nd derivatives.        
c   note the inner dimension of the data is the number of profiles.             
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call spfmax(l,n,x,f,s,p,xp,fp)                                      
c                                                                               
c   input argument list:                                                        
c     l        - integer number of profiles                                     
c     n        - integer number of points in each profile                       
c     x        - real (n) monotonically increasing abscissa values              
c     f        - real (l,n) data values                                         
c     s        - real (l,n) 2nd derivative of f (from subprogram spcoef)        
c                                                                               
c   output argument list:                                                       
c     p        - real (l) point number                                          
c     xp       - real (l) abscissa values of maximum value                      
c     fp       - real (l) maximum data values                                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension x(n),f(l,n),s(l,n),p(l),xp(l),fp(l)                             
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  find maximum given value                                                     
      do i=1,l                                                                  
        p(i)=1                                                                  
        fp(i)=f(i,1)                                                            
      enddo                                                                     
      do k=2,n                                                                  
        do i=1,l                                                                
          if(f(i,k).gt.fp(i)) then                                              
            p(i)=k                                                              
            fp(i)=f(i,k)                                                        
          endif                                                                 
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  determine maximum value of cubic spline                                      
      do i=1,l                                                                  
        k1=nint(p(i))                                                           
        kt=k1+sign(1,n+1-2*k1)                                                  
        dx=x(k1)-x(kt)                                                          
        df=f(i,k1)-f(i,kt)                                                      
        s1=s(i,k1)                                                              
        st=s(i,kt)                                                              
        dp=df/dx+dx*(2*s1+st)/6                                                 
        k2=k1+sign(1.,dp)                                                       
        if(k2.ge.1.and.k2.le.n) then                                            
          x1=x(k1)                                                              
          x2=x(k2)                                                              
          xm=(x2+x1)/2                                                          
          dx=x2-x1                                                              
          f1=f(i,k1)                                                            
          f2=f(i,k2)                                                            
          df=f2-f1                                                              
          s1=s(i,k1)                                                            
          s2=s(i,k2)                                                            
          sm=(s2+s1)/2                                                          
          ds=s2-s1                                                              
          if(ds.ne.0.) then                                                     
            xpa=xm-sm*dx/ds                                                     
            xpb=sqrt((dx**2*(4*sm**2-s1*s2)/(3*ds)-2*df)/ds)                    
            xp(i)=xpa+xpb                                                       
            sp=s1+ds*(xp(i)-x(k1))/dx                                           
            if(sp.gt.0.) xp(i)=xpa-xpb                                          
          elseif(s1.lt.0.) then                                                 
            xp(i)=xm-df/(dx*s1)                                                 
          else                                                                  
            xp(i)=x1                                                            
          endif                                                                 
          dxp=xp(i)-x1                                                          
          p(i)=k1+dxp/dx                                                        
          fp(i)=f1+dxp/dx*(df-(dx-dxp)*(dxp*ds+dx*(2*s1+s2))/6)                 
        else                                                                    
          p(i)=k1                                                               
          xp(i)=x(k1)                                                           
          fp(i)=f(i,k1)                                                         
        endif                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
