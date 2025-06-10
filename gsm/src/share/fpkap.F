      function fpkap(p)                                                         
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$   subprogram  documentation  block                                         
c                                                                               
c subprogram: fpkap        raise surface pressure to the kappa power.           
c   author: phillips         org: w/nmc2x2   date: 29 dec 82                    
c                                                                               
c abstract: raise surface pressure over 100 kpa to the kappa power              
c   using a rational weighted chebyshev approximation.                          
c   the numerator is of order 2 and the denominator is of order 4.              
c   the pressure range is 40-110 kpa and kappa is defined in fpkapx.            
c   the coeffiecients are set by calling gpkap or including bdpkap.             
c   the accuracy of this approximation is almost 8 decimal places.              
c   on the cray, fpkap is over 10 times faster than exact calculation.          
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             standardized kappa,                           
c                                 increased range and accuracy                  
c                                                                               
c usage:  pkap=fpkap(p)                                                         
c                                                                               
c   input argument list:                                                        
c     p        - real surface pressure in kilopascals (cb)                      
c                p should be in the range 40. to 110.                           
c                                                                               
c   output argument list:                                                       
c     fpkap    - real p/100 to the kappa power                                  
c                                                                               
c common blocks:                                                                
c   compkap  - coefficients for function fpkap                                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      common/compkap/ cn0,cn1,cn2,cd0,cd1,cd2,cd3,cd4                           
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      fpkap=(cn0+p*(cn1+p*cn2))/(cd0+p*(cd1+p*(cd2+p*(cd3+p*cd4))))             
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
