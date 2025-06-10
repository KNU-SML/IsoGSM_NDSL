      function ftlcl(t,tdpd)                                                    
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$   subprogram  documentation  block                                         
c                                                                               
c subprogram: ftlcl        compute lcl temperature.                             
c   author: phillips         org: w/nmc2x2   date: 29 dec 82                    
c                                                                               
c abstract: compute temperature at the lifting condensation level               
c   from temperature and dewpoint depression.  the formula used is              
c   a polynomial taken from phillips mstadb routine which empirically           
c   approximates the original exact implicit relationship.                      
c   (this kind of approximation is customary (inman, 1969), but                 
c   the original source for this particular one is not yet known. -mi)          
c   its accuracy is about 0.03 kelvin for a dewpoint depression of 30.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c                                                                               
c usage:  tlcl=ftlcl(t,tdpd)                                                    
c                                                                               
c   input argument list:                                                        
c     t        - real temperature in kelvin                                     
c     tdpd     - real dewpoint depression in kelvin                             
c                                                                               
c   output argument list:                                                       
c     ftlcl    - real temperature at the lcl in kelvin                          
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(clcl1= 0.954442e+0,clcl2= 0.967772e-3,                          
     &          clcl3=-0.710321e-3,clcl4=-0.270742e-5)                          
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      ftlcl=t-tdpd*(clcl1+clcl2*t+tdpd*(clcl3+clcl4*t))                         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
