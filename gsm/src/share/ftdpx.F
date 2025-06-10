      function ftdpx(pv)                                                        
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
cfpp$ expand(ftdp,ftdpxg)                                                       
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftdpx        compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute dewpoint temperature from vapor pressure.           
c   an approximate dewpoint temperature for function ftdpxg                     
c   is obtained using ftdp so gtdp must be already called.                      
c   see documentation for ftdpxg for details.                                   
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   tdp=ftdpx(pv)                                                        
c                                                                               
c   input argument list:                                                        
c     pv       - real vapor pressure in kilopascals (cb)                        
c                                                                               
c   output argument list:                                                       
c     ftdpx    - real dewpoint temperature in kelvin                            
c                                                                               
c subprograms called:                                                           
c   (ftdp)   - inlinable function to compute dewpoint temperature               
c   (ftdpxg) - inlinable function to compute dewpoint temperature               
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      tg=ftdp(pv)                                                               
      ftdpx=ftdpxg(tg,pv)                                                       
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
