      function fpvsq(t)                                                         
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fpvsq        compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute saturation vapor pressure from the temperature.             
c   a quadratic interpolation is done between values in a lookup table          
c   computed in gpvs. see documentation for fpvsx for details.                  
c   input values outside table range are reset to table extrema.                
c   the interpolation accuracy is almost 9 decimal places.                      
c   on the cray, fpvsq is about 3 times faster than exact calculation.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             quadratic interpolation                       
c                                                                               
c usage:   pvs=fpvsq(t)                                                         
c                                                                               
c   input argument list:                                                        
c     t        - real temperature in kelvin                                     
c                                                                               
c   output argument list:                                                       
c     fpvsq    - real saturation vapor pressure in kilopascals (cb)             
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
      dimension tbpvs(nx),tbpvs0(nx)
      common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
      common/compvs/ c1xpvs,c2xpvs,tbpvs
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xj=min(max(c1xpvs+c2xpvs*t,1.),float(nx))
      jx=min(max(nint(xj),2),nx-1)
      dxj=xj-jx
      fj1=tbpvs(jx-1)
      fj2=tbpvs(jx)
      fj3=tbpvs(jx+1)
      fpvsq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
