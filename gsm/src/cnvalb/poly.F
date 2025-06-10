      subroutine poly(n,rad,p)                                                  
c     implicit double precision (a-h,o-z)                                       
      x  = cos(rad)                                                             
      y1 = 1.0e+00                                                              
      y2 = x                                                                    
      do 10 i=2,n                                                               
        g  = x * y2                                                             
        y3 = g - y1 + g - (g-y1) / float(i)                                     
        y1 = y2                                                                 
        y2 = y3                                                                 
  10  continue                                                                  
      p = y3                                                                    
      return                                                                    
      end                                                                       
