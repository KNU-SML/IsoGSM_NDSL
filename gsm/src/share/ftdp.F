      function ftdp(pv)                                                         
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftdp         compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute dewpoint temperature from vapor pressure.                   
c   a linear interpolation is done between values in a lookup table             
c   computed in gtdp. see documentation for ftdpxg for details.                 
c   input values outside table range are reset to table extrema.                
c   the interpolation accuracy is better than 0.0005 kelvin                     
c   for dewpoint temperatures greater than 250 kelvin,                          
c   but decreases to 0.02 kelvin for a dewpoint around 230 kelvin.              
c   on the cray, ftdp is about 75 times faster than exact calculation.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:   tdp=ftdp(pv)                                                         
c                                                                               
c   input argument list:                                                        
c     pv       - real vapor pressure in kilopascals (cb)                        
c                                                                               
c   output argument list:                                                       
c     ftdp     - real dewpoint temperature in kelvin                            
c                                                                               
c common blocks:                                                                
c   comtdp   - scaling parameters and table computed in gtdp.                   
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=5001)                                                        
      dimension tbtdp(nx)                                                       
      common/comtdp/ c1xtdp,c2xtdp,tbtdp                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xj=min(max(c1xtdp+c2xtdp*pv,1.),float(nx))                                
      jx=min(xj,nx-1.)                                                          
      ftdp=tbtdp(jx)+(xj-jx)*(tbtdp(jx+1)-tbtdp(jx))                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
