      function ftheq(t,pk)                                                      
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftheq        compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute equivalent potential temperature at the lcl                 
c   from temperature and pressure over 100 kpa to the kappa power.              
c   a biquadratic interpolation is done between values in a lookup table        
c   computed in gthe. see documentation for fthex for details.                  
c   input values outside table range are reset to table extrema,                
c   except zero is returned for too cold or high lcls.                          
c   the interpolation accuracy is better than 0.0002 kelvin.                    
c   on the cray, ftheq is almost 3 times faster than exact calculation.         
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             quadratic interpolation                       
c                                                                               
c usage:   the=ftheq(pv)                                                        
c                                                                               
c   input argument list:                                                        
c     t        - real lcl temperature in kelvin                                 
c     pk       - real lcl pressure over 100 kpa to the kappa power              
c                                                                               
c   output argument list:                                                       
c     ftheq    - real equivalent potential temperature in kelvin                
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
        jx=min(max(nint(xj),2),nx-1)                                            
        jy=min(max(nint(yj),2),ny-1)                                            
        dxj=xj-jx                                                               
        dyj=yj-jy                                                               
        ft11=tbthe(jx-1,jy-1)                                                   
        ft12=tbthe(jx-1,jy)                                                     
        ft13=tbthe(jx-1,jy+1)                                                   
        ft21=tbthe(jx,jy-1)                                                     
        ft22=tbthe(jx,jy)                                                       
        ft23=tbthe(jx,jy+1)                                                     
        ft31=tbthe(jx+1,jy-1)                                                   
        ft32=tbthe(jx+1,jy)                                                     
        ft33=tbthe(jx+1,jy+1)                                                   
        ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21                  
        ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22                  
        ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23                  
        ftheq=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2                 
      else                                                                      
        ftheq=0.                                                                
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
