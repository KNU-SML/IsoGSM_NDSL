      subroutine glats(lats,gaulat)                                             
c     implicit double precision (a-h,o-z)                                       
c     real  gaulat                                                              
      dimension gaulat(lats)                                                    
c                                                                               
      lat2 = lats / 2                                                           
      eps = 1.0e-12                                                             
      si = 1.0e+00                                                              
      k2 = 2 * lat2                                                             
      k1 = k2 - 1                                                               
      pi = atan(si) * 4.0e+00                                                   
      dradz = pi / 360.0e+00                                                    
      rad = 0.0e+00                                                             
      do 100 k=1,lat2                                                           
        drad = dradz                                                            
  10    call poly(k2,rad,p2)                                                    
  20    p1 = p2                                                                 
        rad = rad + drad                                                        
        call poly(k2,rad,p2)                                                    
        if (sign(si,p1) .eq. sign(si,p2)) goto 20                               
        if (drad .lt. eps) goto 30                                              
        rad = rad - drad                                                        
        drad = drad * 0.25e+00                                                  
        go to 10                                                                
  30    continue                                                                
        gaulat(k) = 0.5e+00 * pi - rad                                          
        gaulat(lats-k+1) = -gaulat(k)                                           
 100  continue                                                                  
      if (mod(lats,2) .eq. 1) gaulat(lat2+1) = 0.0                              
      return                                                                    
      end                                                                       
