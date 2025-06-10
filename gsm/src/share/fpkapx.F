#define CP cp
#define RD rd
      function fpkapx(p)                                                        
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$   subprogram  documentation  block                                         
c                                                                               
c subprogram: fpkapx       raise surface pressure to the kappa power.           
c   author: phillips         org: w/nmc2x2   date: 29 dec 82                    
c                                                                               
c abstract: raise surface pressure over 100 kpa to the kappa power.             
c   kappa is equal to rd/cp where rd and cp are physical constants.             
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   94-12-30  iredell             made into inlinable function                  
c                                                                               
c usage:  pkap=fpkapx(p)                                                        
c                                                                               
c   input argument list:                                                        
c     p        - real surface pressure in kilopascals (cb)                      
c                                                                               
c   output argument list:                                                       
c     fpkapx   - real p/100 to the kappa power                                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(cp=cp_,rd=rd_)                                                  
      parameter(rocp=rd/cp)                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      fpkapx=(p/100.)**rocp                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
