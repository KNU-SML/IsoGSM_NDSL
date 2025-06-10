#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine gaeros(si,sl)                                                  
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c********************************************************************           
c  setup common block 'swaer' for aerosols and rayleigh scattering              
c  optical properties in four uv+vis bands and four nir bands.                  
c   band: 1. 0.225-0.285 (uv)       2. 0.175-0.225;0.285-0.300 (uv)             
c         3. 0.300-0.325 (uv0       4. 0.325-0.690 (par)                        
c         5. 2.27 - 4.0  (nir)      6. 1.22 - 2.27 (nir)                        
c         7. 0.70 - 1.22 (nir)      8. 0.70 - 4.0  (nir)                        
c  ref: wmo report wcp-112 (1986)                                               
c                                                                               
c  the six typical aerosol profiles:                                            
c     1.ubr; 2.cont-1; 3.mar-1; 4.cont-2; 5.mar-2; 6.conv                       
c                                                                               
c     sigref  - ref. sigma level                   n/d   ndm*nae                
c  arrays in the common block:                                                  
c     sig0    - ratio of sigl to domain boundary   n/d   l*nae                  
c     haer    - scale height of aerosols           km    ndm*nae                
c     zaer    - ext. coef. of aerosols             1/km  ndm*nae                
c     hh      - atmospheric scale height           km    l                      
c     dz      - layer thickness                    km    l                      
c     idm     - domain index                       n/d   l*nae                  
c     oaer    - single scattering albedo           n/d   ndm*nae*nbd            
c     gaer    - asymmetry parameter                n/d   ndm*nae*nbd            
c     raer    - ratio of band wavelength to the                                 
c               refference wavelength (0.55 micron)n/d   nbd                    
c     taur    - rayleigh scattering optical depth  n/d   l*nbd                  
c********************************************************************           
c                                                                               
#include <paramodel.h>
#include <rdparm.h>
#include <comswaer.h>
      dimension si(lp1),sl(l),sigln(lp1)
c
c===> ... compute layer distributions of rayleigh scattering                    
c                                                                               
      do 10 n=1,nbd                                                             
      do 10 k=1,l                                                               
        taur(k,n) = tauray(n) * (si(lp1-k)-si(lp1-k+1))                         
  10  continue                                                                  
c                                                                               
c===> ... setup log sigma array (set toa sigma=0.0001)                          
c         rem: si,sl k=1 is sfc; but in radiation k=1 is toa                    
      sigln(1) = alog(1.0e-4)                                                  
      do 20 k=1,l                                                               
        sigln(k+1) = alog(si(lp1-k))                                            
  20  continue                                                                  
      do 30 k=1,l                                                               
        hh(k) = 6.05e0 + 2.5e0 * sl(lp1-k)                                    
        dz(k) = hh(k) * (sigln(k+1)-sigln(k))                                   
  30  continue                                                                  
c                                                                               
      do 60 iaer=1,nae                                                          
        sbund= si(1)                                                            
        idom = 1                                                                
        do 50 k=l,1,-1                                                          
          if (si(lp1-k+1) .lt. sigref(idom,iaer)) then                          
            idom = idom + 1                                                     
            sbund= si(lp1-k)                                                    
          end if                                                                
          sig0(k,iaer) = sl(lp1-k) / sbund                                      
          idm (k,iaer) = idom                                                   
  50    continue                                                                
  60  continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
