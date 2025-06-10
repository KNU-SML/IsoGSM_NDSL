      function fpvs(t)                                                          
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fpvs         compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute saturation vapor pressure from the temperature.             
c   a linear interpolation is done between values in a lookup table             
c   computed in gpvs. see documentation for fpvsx for details.                  
c   input values outside table range are reset to table extrema.                
c   the interpolation accuracy is almost 6 decimal places.                      
c   on the cray, fpvs is about 4 times faster than exact calculation.           
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             expand table                                  
c   96-02-19  hong                ice effect
c                                                                               
c usage:   pvs=fpvs(t)                                                          
c                                                                               
c   input argument list:                                                        
c     t        - real temperature in kelvin                                     
c                                                                               
c   output argument list:                                                       
c     fpvs     - real saturation vapor pressure in kilopascals (cb)             
c                                                                               
c common blocks:                                                                
c   compvs   - scaling parameters and table computed in gpvs.                   
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=7501)                                                        
      dimension tbpvs(nx)                                                       
      common/compvs/ c1xpvs,c2xpvs,tbpvs                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xj=min(max(c1xpvs+c2xpvs*t,1.),float(nx))                                 
      jx=min(xj,nx-1.)                                                          
      fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
