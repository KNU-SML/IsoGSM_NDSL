#define G g
#define RD rd
      subroutine hydrosp(im,km,si,sl,zs,tv,z)                              
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    hydrosp     calculate geopotential heights                     
c                                                                               
c abstract: calculates spectral coeffs of geopotential heights 
c on the sigma full levels as a function of orography, virtual temperature           
c     dz=rd/g*tv*dlnp                                                           
c                                                                               
c program history log:                                                          
c                                                                               
c usage:    call hydro(im,km,si,sl,zs,tv,q,z)                              
c                                                                               
c   input argument list:                                                        
c     im       - integer number of coefficients
c     km       - integer number of levels                                       
c     si       - real (km+1) sigma interface values                             
c     sl       - real (km) sigma values                                         
c     zs       - real (im) orography in m
c     tv       - real (im,km) vertual temperature in k                                  
c                                                                               
c   output argument list:                                                       
c     z        - real (im,km) heights on the full levels in m                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
#include <paramodel.h>
c$$$                                                                            
      dimension si(km+1),sl(km),zs(im),tv(im,km)
      dimension z(im,km)
      parameter(rog_=rd_/g_)                                        
c
      dimension zi(im,km)
c
      do i=1,im                                                                 
        zi(i,1)=zs(i)                                                           
      enddo                                                                     
c
      do k=1,km-1                                                               
        ca=rog_*log(si(k)/sl(k))                                                 
        cb=rog_*log(sl(k)/si(k+1))                                               
        do i=1,im                                                               
          z(i,k)=zi(i,k)+ca*tv(i,k)
          zi(i,k+1)=z(i,k)+cb*tv(i,k)
        enddo                                                                   
      enddo                                                                     
c
      ca=rog_*log(si(km)/sl(km))                                                 
      do i=1,im                                                                 
        z(i,km)=zi(i,km)+ca*tv(i,k)
      enddo                                                                     
c
      return                                                                    
      end                                                                       
