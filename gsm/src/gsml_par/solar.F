#define PI pi
#undef sun
      subroutine solar(jd,fjd,r,dlt,alp,slag,sdec,cdec)                         
c>yh  subroutine solar(jd,fjd,r,dlt,alp,slag)                                   
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    solar       astronomical(solar) data - sw radiation            
c   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07                     
c                                                                               
c abstract: this code written at gfdl ....                                      
c   computes radius vector,declination and right ascension of sun,              
c   and equation of time. same as subroutine -solmrf-, but                      
c   with hour angle,fractional daylight,and                                     
c   mean zenith angle calculations removed (subroutine -zenith-                 
c   calculates these for each point rather than each latitude).                 
c   --this code is to be used for other fcst models or for the mrf              
c   --model if 'instantaneous' sw calculations desired.                         
c                                                                               
c program history log:                                                          
c   77-07-21  robert white,gfdl.                                                
c   89-07-07  kenneth campana-moved the hour angle calculations to              
c                              subroutine -zenith-                              
c                                                                               
c usage:    call solmrf(jd,fjd,r,dlt,alp,slag)                                  
c   input argument list:                                                        
c     jd       - julian day for current fcst hour.                              
c     fjd      - fraction of the julian day.                                    
c   output argument list:                                                       
c     r        - radius vector of the sun.                                      
c     dlt      - declination of the sun (radians).                              
c     alp      - right ascension of the sun.                                    
c     slag     - equation of time (radians).                                    
c                                                                               
c   output files:                                                               
c     output   - print file.                                                    
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter (pi=pi_,tpi=2.0*pi,hpi=0.5*pi,rad=180.0/pi)                     
c                                                                               
c                                                                               
c    *******************************************************************        
c    *                            s o l a r                            *        
c... *  patterned after original gfdl code---                          *        
c... *     but no calculation of latitude mean cos solar zenith angle..*        
c... *     zenith angle calculations done in subr zenith in this case..*        
c... *  hr angle,mean cosz,and mean tauda calc removed--k.a.c. mar 89  *        
c    *  updates by hualu pan to limit iterations in newton method and  *        
c    *  also ccr reduced from(1.3e-7)--both to avoid nonconvergence in *        
c    *  nmc s half precision version of gfdl s code   ----  fall 1988  *        
c    *******************************************************************        
c                                                                               
c.....solar computes radius vector, declination and right ascension of          
c.....sun, equation of time                                                     
c                                                                               
                                   d a t a                                      
     1   cyear/365.25/,      ccr/1.3e-6/                                       
c                                                                               
c.....tpp = days between epoch and perihelion passage of 1900                   
c.....svt6 = days between perihelion passage and march equinox of 1900          
c.....jdor = jd of epoch which is january 0, 1900 at 12 hours ut                
c                                                                               
                                   d a t a                                      
     1   tpp/1.55/,          svt6/78.035/,       jdor/2415020/                  
c                                                                               
c    *******************************************************************        
c                                                                               
      dat=float(jd-jdor)-tpp+fjd                                                
c    computes time in julian centuries after epoch                              
      t=float(jd-jdor)/36525.e0                                                
c    computes length of anomalistic and tropical years (minus 365 days)         
      year=.25964134e0+.304e-5*t                                              
      tyear=.24219879e0-.614e-5*t                                             
c    computes orbit eccentricity and angle of earths inclination from t        
      ec=.01675104e0-(.418e-4+.126e-6*t)*t                                   
      angin=23.452294e0-(.0130125e0+.164e-5*t)*t                             
      ador=jdor                                                                 
      jdoe=ador+(svt6*cyear)/(year-tyear)                                       
c    deleqn=updated svt6 for current date                                       
      deleqn=float(jdoe-jd)*(year-tyear)/cyear                                  
      year=year+365.e0                                                         
      sni=sin(angin/rad)                                                        
      tini=1.e0/tan(angin/rad)                                                 
      er=sqrt((1.e0+ec)/(1.e0-ec))                                            
      qq=deleqn*tpi/year                                                        
c    determine true anomaly at equinox                                          
      e=1.e0                                                                   
      iter = 0                                                                  
 32   ep=e-(e-ec*sin(e)-qq)/(1.e0-ec*cos(e))                                   
      cd=abs(e-ep)                                                              
      e=ep                                                                      
      iter = iter + 1                                                           
#ifndef NOPRINT
      if(iter.gt.10) then                                                       
        write(6,*) ' iteration count for loop 32 =', iter                       
        write(6,*) ' e, ep, cd =', e, ep, cd                                    
      endif                                                                     
#endif
      if(iter.gt.10) goto 1032                                                  
      if(cd.gt.ccr) go to 32                                                    
 1032 continue                                                                  
      he=.5e0*e                                                                
      eq=2.e0*atan(er*tan(he))                                                 
c    date=days since last perihelion passage                                    
      date = mod(dat,year)                                                      
c    solve orbit equations by newtons method                                   
      em=tpi*date/year                                                          
      e=1.e0                                                                   
      iter = 0                                                                  
 31   ep=e-(e-ec*sin(e)-em)/(1.e0-ec*cos(e))                                   
      cr=abs(e-ep)                                                              
      e=ep                                                                      
      iter = iter + 1                                                           
#ifndef NOPRINT
      if(iter.gt.10) then                                                       
        write(6,*) ' iteration count for loop 31 =', iter                       
      endif                                                                     
#endif
      if(iter.gt.10) goto 1031                                                  
      if(cr.gt.ccr) go to 31                                                    
 1031 continue                                                                  
      r=1.e0-ec*cos(e)                                                         
      he=.5e0*e                                                                
      w=2.e0*atan(er*tan(he))                                                  
c>yh  sind=sni*sin(w-eq)                                                        
c>yh  dlt=asin(sind)                                                            
      sdec=sni*sin(w-eq)                                                        
      cdec=sqrt(1.e0 - sdec*sdec)                                              
      dlt=asin(sdec)                                                            
      alp=asin(tan(dlt)*tini)                                                   
      tst=cos(w-eq)                                                             
      if(tst.lt.0.e0) alp=pi-alp                                               
      if(alp.lt.0.e0) alp=alp+tpi                                              
      sun=tpi*(date-deleqn)/year                                                
      if(sun.lt.0.e0) sun=sun+tpi                                              
      slag=sun-alp-.03255e0                                                    
      return                                                                    
      end                                                                       
