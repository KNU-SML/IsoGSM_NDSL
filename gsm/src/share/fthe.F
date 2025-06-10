      function fthe(t,pk)                                                       
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fthe         compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute equivalent potential temperature at the lcl                 
c   from temperature and pressure over 100 kpa to the kappa power.              
c   a bilinear interpolation is done between values in a lookup table           
c   computed in gthe. see documentation for fthex for details.                  
c   input values outside table range are reset to table extrema,                
c   except zero is returned for too cold or high lcls.                          
c   the interpolation accuracy is better than 0.01 kelvin.                      
c   on the cray, fthe is almost 6 times faster than exact calculation.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:   the=fthe(pv)                                                         
c                                                                               
c   input argument list:                                                        
c     t        - real lcl temperature in kelvin                                 
c     pk       - real lcl pressure over 100 kpa to the kappa power              
c                                                                               
c   output argument list:                                                       
c     fthe     - real equivalent potential temperature in kelvin                
c                                                                               
c common blocks:                                                                
c   comthe   - scaling parameters and table computed in gthe.                   
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=241,ny=151)                                                  
      dimension tbthe(nx,ny)                                                    
      common/comthe/ c1xthe,c2xthe,c1ythe,c2ythe,tbthe                          
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xj=min(c1xthe+c2xthe*t,float(nx))                                         
      yj=min(c1ythe+c2ythe*pk,float(ny))                                        
      if(xj.ge.1..and.yj.ge.1.) then                                            
        jx=min(xj,nx-1.)                                                        
        jy=min(yj,ny-1.)                                                        
        ftx1=tbthe(jx,jy)+(xj-jx)*(tbthe(jx+1,jy)-tbthe(jx,jy))                 
        ftx2=tbthe(jx,jy+1)+(xj-jx)*(tbthe(jx+1,jy+1)-tbthe(jx,jy+1))           
        fthe=ftx1+(yj-jy)*(ftx2-ftx1)                                           
      else                                                                      
        fthe=0.                                                                 
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
