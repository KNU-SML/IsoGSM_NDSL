      subroutine poly(n,rad,p)                                                  
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    poly        evaluates legendre polynomial.                     
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-01                     
c                                                                               
c abstract: evaluates the unnormalized legendre polynomial                      
c   of specified degree at a given colatitude using a standard                  
c   recursion formula.  real arithmetic is used.                                
c                                                                               
c program history log:                                                          
c   88-04-01  joseph sela                                                       
c                                                                               
c usage:    call poly (n, rad, p)                                               
c   input argument list:                                                        
c     n        - degree of legendre polynomial.                                 
c     rad      - real colatitude in radians.                                    
c                                                                               
c   output argument list:                                                       
c     p        - real value of legendre polynomial.                             
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      x = cos(rad)                                                              
      y1 = 1.0                                                                  
      y2=x                                                                      
      do 1 i=2,n                                                                
      g=x*y2                                                                    
      y3=g-y1+g-(g-y1)/float(i)                                                 
      y1=y2                                                                     
      y2=y3                                                                     
1     continue                                                                  
      p=y3                                                                      
      return                                                                    
      end                                                                       
