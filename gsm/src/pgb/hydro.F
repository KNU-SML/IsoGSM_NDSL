      subroutine hydro(im,ix,km,si,sl,zs,t,q,z,zi)                              
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    hydro       calculate geopotential heights                     
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: calculates geopotential heights on both the sigma interfaces        
c   and the sigma full levels as a function of orography, temperature           
c   and moisture.  virtual temperature is calculated from temperature           
c   and moisture and the hydrostatic equation is integrated                     
c     dz=rd/g*tv*dlnp                                                           
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call hydro(im,ix,km,si,sl,zs,t,q,z,zi)                              
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of levels                                       
c     si       - real (km+1) sigma interface values                             
c     sl       - real (km) sigma values                                         
c     zs       - real (im) orography is m                                       
c     t        - real (ix,km) temperature in k                                  
c     q        - real (ix,km) specific humidity in kg/kg                        
c                                                                               
c   output argument list:                                                       
c     z        - real (ix,km) heights on the full levels in m                   
c     zi       - real (ix,km) heights on the interfaces in m                    
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension si(ix,km+1),sl(ix,km),zs(im),t(ix,km),q(ix,km)                        
      dimension z(ix,km),zi(ix,km)                                              
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 ,rv= 4.6150e+2 )                    
      parameter(rog=rd/g,fvirt=rv/rd-1.)                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do i=1,im                                                                 
        zi(i,1)=zs(i)                                                           
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do k=1,km-1                                                               
        do i=1,im                                                               
          ca=rog*log(si(i,k)/sl(i,k))                                                 
          cb=rog*log(sl(i,k)/si(i,k+1))                                               
          tv=t(i,k)*(1.+fvirt*q(i,k))                                           
          z(i,k)=zi(i,k)+ca*tv                                                  
          zi(i,k+1)=z(i,k)+cb*tv                                                
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do i=1,im                                                                 
        ca=rog*log(si(i,km)/sl(i,km))                                                 
        tv=t(i,km)*(1.+fvirt*q(i,km))                                           
        z(i,km)=zi(i,km)+ca*tv                                                  
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
