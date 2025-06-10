      function ftmax(the,pk,qma)                                                
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
cfpp$ expand(ftma,ftmaxg)                                                       
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftmax        compute moist adiabat temperature                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute temperature and humidity of a parcel                
c   lifted up a moist adiabat from equivalent potential temperature             
c   at the lcl and pressure over 100 kpa to the kappa power.                    
c   an approximate parcel temperature for function ftmaxg                       
c   is obtained using ftma so gtma must be already called.                      
c   see documentation for ftmaxg for details.                                   
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   tma=ftmax(the,pk,qma)                                                
c                                                                               
c   input argument list:                                                        
c     the      - real equivalent potential temperature in kelvin                
c     pk       - real pressure over 100 kpa to the kappa power                  
c                                                                               
c   output argument list:                                                       
c     ftmax    - real parcel temperature in kelvin                              
c     qma      - real parcel specific humidity in kg/kg                         
c                                                                               
c subprograms called:                                                           
c   (ftma)   - inlinable function to compute parcel temperature                 
c   (ftmaxg) - inlinable function to compute parcel temperature                 
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      tg=ftma(the,pk,qg)                                                        
      ftmax=ftmaxg(tg,the,pk,qma)                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
