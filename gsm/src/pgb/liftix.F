#define IO io
#define LEVS levs
      subroutine liftix(im,ix,km,sl,ps,t,q,tm,tli)                              
c$$$  subprogram documentation block                                            
c                                                                               
c  subprogram:    liftix      compute best lifted index from sigma              
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: computes the best lifted index from profiles in sigma.              
c   the best lifted index is here computed by finding the parcel                
c   below sigma 0.80 with the warmest equivalent potential temperature,         
c   then raising it to 500 mb and subtracting its parcel temperature            
c   from the environment temperature.                                           
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c   94-04-28  iredell   fixed parameters                                        
c                                                                               
c usage:    call liftix(im,ix,km,sl,ps,t,q,tm,tli)                              
c                                                                               
c   input argument list:                                                        
c     im       - integer number of points                                       
c     ix       - integer first dimension of upper air data                      
c     km       - integer number of sigma levels                                 
c     sl       - real (km) sigma values                                         
c     ps       - real (im) surface pressure in kpa                              
c     t        - real (ix,km) temperature in k                                  
c     q        - real (ix,km) specific humidity in kg/kg                        
c     tm       - real (im) 500 mb temperature in k                              
c                                                                               
c   output argument list:                                                       
c     tli      - real (ix,km) best lifted index in k                            
c                                                                               
c subprograms called:                                                           
c   (fpkap)   - function to compute pressure to the kappa                       
c   (ftdp)    - function to compute dewpoint temperature                        
c   (ftlcl)   - function to compute lifting condensation level                  
c   (fthe)    - function to compute equivalent potential temperature            
c   (ftma)    - function to compute moist adiabat temperature                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#ifdef CRAY_THREAD
cfpp$ expand(fpkap,ftdp,ftlcl,fthe,ftma)                                        
#endif
#include <paramodel.h>
      dimension sl(km),ps(im),t(ix,km),q(ix,km),tm(im),tli(im)                  
      parameter(cp= 1.0046e+3 ,rd= 2.8705e+2 ,rv= 4.6150e+2 )                   
      parameter(rk=rd/cp,eps=rd/rv,epsm1=rd/rv-1.)                              
      parameter(slift=0.80,plift=50.)                                           
      dimension slk(levs_),psk(2*io_),slkma(2*io_),thema(2*io_)                 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  initialize variables                                                         
      pliftk=(plift/100.)**rk                                                   
      do k=1,km                                                                 
        slk(k)=sl(k)**rk                                                        
      enddo                                                                     
      do i=1,im                                                                 
        psk(i)=fpkap(ps(i))                                                     
        slkma(i)=0.                                                             
        thema(i)=0.                                                             
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  select the warmest equivalent potential temperature                          
c  between the surface and sigma slift1                                         
      k=1                                                                       
      dowhile(sl(k).gt.slift)                                                   
        do i=1,im                                                               
          p=sl(k)*ps(i)                                                         
          pv=p*q(i,k)/(eps-epsm1*q(i,k))                                        
          tdpd=max(t(i,k)-ftdp(pv),0.)                                          
          tlcl=ftlcl(t(i,k),tdpd)                                               
          slklcl=slk(k)*tlcl/t(i,k)                                             
          thelcl=fthe(tlcl,slklcl*psk(i))                                       
          if(thelcl.gt.thema(i)) then                                           
            slkma(i)=slklcl                                                     
            thema(i)=thelcl                                                     
          endif                                                                 
        enddo                                                                   
        k=k+1                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  lift the parcel to 500 mb along a dry adiabat below the lcl                  
c  or along a moist adiabat above the lcl.                                      
c  the lifted index is the environment minus parcel temperature.                
      do i=1,im                                                                 
        if(ps(i).gt.plift.and.thema(i).gt.0.) then                              
          slkp=pliftk/psk(i)                                                    
          slkc=min(slkp,slkma(i))                                               
          tlift=slkp/slkc*ftma(thema(i),slkc*psk(i),qma)                        
          tli(i)=tm(i)-tlift                                                    
        else                                                                    
          tli(i)=0.                                                             
        endif                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
