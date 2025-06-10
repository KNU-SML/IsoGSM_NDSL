      subroutine getrh(im,ix,km,sl,ps,q,t,r)                                    
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
cfpp$ expand(fpvs,fpvs0)                                                              
#endif
#include <paramodel.h>
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    getrh       calculate relative humidity                        
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: calculates relative humidity as a function of pressure,             
c   specific humidity and moisture.  saturation specific humidity               
c   is calculated from saturation vapor pressure which is returned              
c   from a lookup table routine fpvs.                                           
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call getrh(im,ix,km,sl,ps,q,t,r)                                    
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of levels                                       
c     sl       - real (km) sigma values                                         
c     ps       - real (im) surface pressure in kpa                              
c     q        - real (ix,km) specific humidity in kg/kg                        
c     t        - real (ix,km) temperature in k                                  
c                                                                               
c   output argument list:                                                       
c     r        - real (ix,km) relative humidity in percent                      
c                                                                               
c subprograms called:                                                           
c   (fpvs)   - function to compute saturation vapor pressure                    
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension sl(ix,km),ps(im),q(ix,km),t(ix,km),r(ix,km)                        
      parameter(rd= 2.8705e+2 ,rv= 4.6150e+2 ,eps=rd/rv,epsm1=rd/rv-1.)         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do k=1,km                                                                 
        do i=1,im                                                               
c !!! for comparison with R-2, rh is calculated without ICE
c#ifdef ICE
c          es=fpvs(t(i,k))                                                       
c#else
          es=fpvs0(t(i,k))                                                       
c#endif
          qs=eps*es/(sl(i,k)*ps(i)+epsm1*es)                                      
          r(i,k)=min(max(q(i,k)/qs*100.,0.),100.)                               
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
