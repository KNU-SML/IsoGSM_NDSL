#define LEVS levs
      subroutine spcoef(l,n,x,f,s)                                              
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    spcoef      compute 2nd derivatives for cubic splines         
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: compute the second derivatives of cubic spline profiles             
c   in preparation for cubic spline interpolations.                             
c   cubic splines are piecewise cubic polynomials fitting the data              
c   with continuous first and second derivatives at interior points             
c   and second derivatives set to zero at and beyond the end points.            
c   the computations are done by marching up then down the profiles.            
c   note the inner dimension of the data is the number of profiles.             
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call spcoef(l,n,x,f,s)                                              
c                                                                               
c   input argument list:                                                        
c     l        - integer number of profiles                                     
c     n        - integer number of points in each profile                       
c     x        - real (n) monotonically increasing abscissa values              
c     f        - real (l,n) data values                                         
c                                                                               
c   output argument list:                                                       
c     s        - real (l,n) 2nd derivative of f with respect to x               
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension x(l,n),f(l,n),s(l,n)                                              
      dimension rho(levs_-1)                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  initialize end points                                                        
      rho(1)=0.                                                                 
      do i=1,l                                                                  
        s(i,1)=0.                                                               
        s(i,n)=0.                                                               
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  march up the profiles                                                        
      do k=2,n-1                                                                
        do i=1,l                                                                
          hm1=x(i,k)-x(i,k-1)                                                         
          rh=1./(x(i,k+1)-x(i,k))                                                     
          rho(k)=-1./(hm1*(rho(k-1)+2.)*rh+2.)                                    
          d=6.*((f(i,k+1)-f(i,k))*rh-(f(i,k)-f(i,k-1))/hm1)*rh                  
          s(i,k)=(hm1*s(i,k-1)*rh-d)*rho(k)                                     
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  march down the profiles                                                      
      do k=n-1,2,-1                                                             
        do i=1,l                                                                
          s(i,k)=rho(k)*s(i,k+1)+s(i,k)                                         
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
