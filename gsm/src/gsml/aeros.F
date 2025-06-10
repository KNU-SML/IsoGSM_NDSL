#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine aeros(ipts,ib,kaer,paer,tau,ssa,asy,taurs)                          
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c********************************************************************           
c  compute aerosols optical properties of six typical profiles                  
c  in four uv+vis bands and four nir bands.                                     
c   band: 1. 0.225-0.285 (uv)       2. 0.175-0.225;0.285-0.300 (uv)             
c         3. 0.300-0.325 (uv)       4. 0.325-0.690 (par)                        
c         5. 2.27 - 4.0  (nir)      6. 1.22 - 2.27 (nir)                        
c         7. 0.70 - 1.22 (nir)      8. 0.70 - 4.0  (nir)                        
c  ref: wmo report wcp-112 (1986)                                               
c                                                                               
c  formulations:                                                                
c     dtau(k) = zk(k) * dz(k)                                                   
c     zk(k)   = zk(0) * (prss(k)/press(0))**(h/hd)   for exp type               
c             = const                                for  others                
c     z(k)    = -hh * (ln(press(k))-ln(press(k+1)))                             
c     hh      - atmospheric scale height, a fuction of press                    
c     hd      - aerosols scale height                                           
c   where zk is ext. coeff.; z is height; the vertical indeces are              
c     k=1 at surface for input sigma levels and k=1 at top for all              
c     other quantities.                                                         
c                                                                               
c  input parameters:                                                            
c     ib,kaer,paer                                                              
c                                                                               
c  output parameters:                                                           
c     tau  - optical depth                         n/d                          
c     ssa  - single scattering albedo              n/d                          
c     asy  - asymmetry parameter                   n/d                          
c     taurs- rayleigh scattering optical depth     n/d                          
c********************************************************************           
c                                                                               
#include <paramodel.h>
#include <rdparm.h>
#include <comswaer.h>
c --- input                                                                     
      dimension  paer(imbx,nae)                                                 
c --- output                                                                    
      dimension                                     
     1  tau(imbx,l),  ssa(imbx,l),  asy(imbx,l),  taurs(l)                      
      logical    laer                                                           
      data  crt1,crt2 / 30.0, 0.03333 /                                         
c                                                                               
c===> ... layer optical depth due to rayleigh scattering                        
c                                                                               
      do 20 k=1,l                                                               
        taurs(k) = taur(k,ib)                                                   
  20  continue                                                                  
c                                                                               
      do 30 k=1,l                                                               
      do 30 i=1,ipts                                                            
        ssa(i,k) = 0.0e0                                                       
        asy(i,k) = 0.0e0                                                       
        tau(i,k) = 0.0e0                                                       
  30  continue                                                                  
c                                                                               
      if (kaer .lt. 1) return                                                   
c                                                                               
      do 100 iaer=1,nae                                                         
c                                                                               
        laer = .false.                                                          
        do 40 i=1,ipts                                                          
          laer = laer .or. paer(i,iaer).gt.0.0e0                               
  40    continue                                                                
        if (.not. laer) go to 100                                               
c                                                                               
c===> ... find aerosol optical depth, single scattering albedo                  
c         and asymmetry factor                                                  
        do 80 k=1,l                                                             
          kk = idm(k,iaer)                                                      
          hd = haer(kk,iaer)                                                    
          zk = zaer(kk,iaer)                                                    
          if (hd .gt. 0.0e0) then                                              
            tau0 = zk * (sig0(k,iaer)**(hh(k)/hd)) * dz(k)                      
          else                                                                  
            tau0 = (zk - hd*hh(k)*alog(sig0(k,iaer))) * dz(k)                   
          end if                                                                
          do 60 i=1,ipts                                                        
            tau(i,k) = tau(i,k) + paer(i,iaer)*raer(ib)*tau0                    
            ssa(i,k) = ssa(i,k) + paer(i,iaer)*oaer(kk,iaer,ib)                 
            asy(i,k) = asy(i,k) + paer(i,iaer)*gaer(kk,iaer,ib)                 
  60      continue                                                              
  80    continue                                                                
c                                                                               
 100  continue                                                                  
c===> ... smooth profile at domain boundaries                                   
      do 120 k=2,l                                                              
      do 120 i=1,ipts                                                           
        ratio = 1.0e0                                                          
        if (tau(i,k) .gt. 0.0e0) ratio = tau(i,k-1) / tau(i,k)                 
        tt = tau(i,k) + tau(i,k-1)                                              
        if (ratio .gt. crt1) then                                               
          tau(i,k) = 0.2e0 * tt                                                
          tau(i,k-1) = tt - tau(i,k)                                            
        else if (ratio .lt. crt2) then                                          
          tau(i,k) = 0.8e0 * tt                                                
          tau(i,k-1) = tt - tau(i,k)                                            
        end if                                                                  
 120  continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
