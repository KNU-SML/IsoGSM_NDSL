#define JCAP1 jcap1
#define JCAP jcap
#define LLN22P lln22p
#define LLNP llnp
#define LNT22 lnt22
#define RERTH rerth
      subroutine deldfsp(q,qlap,llstr,llens,lwvdef)
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    deldfsp     compute modified laplacian in spectral space       
c   prgmmr: kanamitsu         org: w/np51     date: 97-03-10                    
c                                                                               
c abstract: computes the modified horizontal laplacian in spectral space        
c           for horizontal diffusion calculation of temp and hum on             
c           grid space.                                                         
c           if n is the total wavenumber, modified laplacian is expressed by:   
c             -(n-n0)*(n-n0+1)/a**2*spc(n,m), n-n0>0                            
c           where n0=0.55*jcap                                                  
c                                                                               
c program history log:                                                          
c   97-03-10  masao kanamitsu                                                   
c                                                                               
c usage:    call deldfsp(q,rq,llstr,llens,lwvdef)                                                  
c                                                                               
c   input argument list:                                                        
c     q        - real (lnt22) ln(ps) spectral coeffs                           
c                                                                               
c   output argument list:                                                       
c     qlap     - real (lnt22) qusi-laplacian of ln(ps) spectral coeffs         
c                                                                               
c attributes:                                                                   
c   language: fortran 77                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
#ifdef MP
#define LLN22S lln22p_
#else
#define LLN22S lnt22_
#endif
      real q(LLN22S)
      real qlap(LLN22S)
      dimension lwvdef(jcap1_)
c                                                                               
      parameter(difcof=3.e15,lefres=80,jdel=2)                                  
c                                                                               
      np=jcap_                                                                  
      n0=0.55*jcap_                                                             
      rtnp=difcof/(rerth_**4)*float(lefres*(lefres+1))**2                       
      jdelh=jdel/2                                                              
      npd=max(np-n0,0)                                                          
      dn1=2.*rtnp/float(npd*(npd+1))**jdelh                                     
c                                                                               
      i=0                                                                       
      do ll=1,llens
        nm=lwvdef(llstr+ll)
        do mm=0,jcap_-nm                                                        
          nd=max(nm+mm-n0,0)                                                    
          dn=dn1*float(nd*(nd+1))**jdelh                                        
          qlap(i+1)=q(i+1)*dn                                                   
          qlap(i+2)=q(i+2)*dn                                                   
          i=i+2                                                                 
        enddo                                                                   
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
