#include <define.h>
      subroutine sunprm(km,ko,po,nsun,                                          
     &                  ipusun,itlsun,ip1sun,ip2sun,kslp)                       
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    sunprm      set parameters for sundry fields                   
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: sets parameters for the sundry fields.                              
c   parameters returned are parameter indicator, level type indicator,          
c   two level numbers and decimal scaling all required for the pds              
c   section of the grib1 message as well relevant sigma layer numbers           
c   for the three lower level relative humidity fields.                         
c   the current nsun=22 sundry fields are:                                      
c     1) surface pressure                                                       
c     2) precipitable water                                                     
c     3) surface orography                                                      
c     4) sea level pressure                                                     
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call sunprm(km,ko,po,nsun,                                          
c    &                  ipusun,itlsun,ip1sun,ip2sun,kslp)                       
c                                                                               
c   input argument list:                                                        
c     km       - integer number of levels                                       
c     ko       - integer number of pressure levels                              
c     po       - real (ko) pressure in millibars                                
c                                                                               
c   output argument list:                                                       
c     ipusun   - integer (nsun) parameter indicators                            
c     itlsun   - integer (nsun) level type indicators                           
c     ip1sun   - integer (nsun) first level numbers                             
c     ip2sun   - integer (nsun) second level numbers                            
c     kslp     - integer (2) relevant pressure levels for slp                   
c                                                                               
c subprograms called:                                                           
c   isrchflex - find first value in an array le target value                     
c   isrcheqx  - find first value in an array equal to target value               
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension po(ko)                                                          
      dimension ipusun(nsun)                                                    
      dimension itlsun(nsun),ip1sun(nsun),ip2sun(nsun)                          
c                                                                               
      dimension kslp(2)                                                         
      dimension pslp(2)                                                         
      data pslp/1000.,500./                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      kslp(1)=mod(isrcheqx(ko,po,1,pslp(1)),ko+1)                               
      kslp(2)=mod(isrcheqx(ko,po,1,pslp(2)),ko+1)                               
      do n=1,nsun                                                               
      ip1sun(n)=0                                                               
      ip2sun(n)=0                                                               
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
