      function fpvs(t)                                                          
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fpvs         compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute saturation vapor pressure from the temperature.             
c   a linear interpolation is done between values in a lookup table             
c   computed in gpvs. see documentation for gpvs for details.                   
c   input values outside table range are reset to table extrema.                
c   this function can be expanded inline in calling routine.                    
c                                                                               
c usage:   pvs=fpvs(t)                                                          
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
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
c$$$                                                                            
      parameter(nx=1501)                                                        
      dimension tbpvs(nx)                                                       
      common/compvs/ c1xpvs,c2xpvs,anxpvs,tbpvs                                 
c-----------------------------------------------------------------------        
      xj=min(max(c1xpvs+c2xpvs*t,1.),anxpvs)                                    
      jx=xj                                                                     
      fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))                            
      return                                                                    
      end                                                                       
