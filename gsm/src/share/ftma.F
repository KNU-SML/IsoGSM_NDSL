      function ftma(the,pk,qma)                                                 
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftma         compute moist adiabat temperature                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute temperature and specific humidity of a parcel               
c   lifted up a moist adiabat from equivalent potential temperature             
c   at the lcl and pressure over 100 kpa to the kappa power.                    
c   bilinear interpolations are done between values in a lookup table           
c   computed in gtma. see documentation for ftmaxg for details.                 
c   input values outside table range are reset to table extrema.                
c   the interpolation accuracy is better than 0.01 kelvin                       
c   and 5.e-6 kg/kg for temperature and humidity, respectively.                 
c   on the cray, ftma is about 35 times faster than exact calculation.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:   tma=ftma(the,pk,qma)                                                 
c                                                                               
c   input argument list:                                                        
c     the      - real equivalent potential temperature in kelvin                
c     pk       - real pressure over 100 kpa to the kappa power                  
c                                                                               
c   output argument list:                                                       
c     ftma     - real parcel temperature in kelvin                              
c     qma      - real parcel specific humidity in kg/kg                         
c                                                                               
c common blocks:                                                                
c   comma    - scaling parameters and table computed in gtma.                   
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=151,ny=121)                                                  
      dimension tbtma(nx,ny),tbqma(nx,ny)                                       
      common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma                         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xj=min(max(c1xma+c2xma*the,1.),float(nx))                                 
      yj=min(max(c1yma+c2yma*pk,1.),float(ny))                                  
      jx=min(xj,nx-1.)                                                          
      jy=min(yj,ny-1.)                                                          
      ftx1=tbtma(jx,jy)+(xj-jx)*(tbtma(jx+1,jy)-tbtma(jx,jy))                   
      ftx2=tbtma(jx,jy+1)+(xj-jx)*(tbtma(jx+1,jy+1)-tbtma(jx,jy+1))             
      ftma=ftx1+(yj-jy)*(ftx2-ftx1)                                             
      qx1=tbqma(jx,jy)+(xj-jx)*(tbqma(jx+1,jy)-tbqma(jx,jy))                    
      qx2=tbqma(jx,jy+1)+(xj-jx)*(tbqma(jx+1,jy+1)-tbqma(jx,jy+1))              
      qma=qx1+(yj-jy)*(qx2-qx1)                                                 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
