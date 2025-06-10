      subroutine compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)                          
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    compjd      computes julian day and fraction                   
c   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07                     
c                                                                               
c abstract: this code written at gfdl ....                                      
c   computes julian day and fraction                                            
c   from year, month, day and time ut...accurate only between                   
c   march 1, 1900 and february 28, 2100.. based on julian calendar              
c   corrected to correspond to gregorian calendar during this period.           
c                                                                               
c program history log:                                                          
c   77-05-06  ray orzol,gfdl                                                    
c                                                                               
c usage:    call compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)                          
c   input argument list:                                                        
c     jyr      - year (4 digits)-intial fcst time.                              
c     jmnth    - month-initial fcst time.                                       
c     jday     - day-initial fcst time.                                         
c     jhr      - z-time of initial fcst time.                                   
c     jmn      - minutes (zero passed from calling program).                    
c   output argument list:                                                       
c     jd       - julian day.                                                    
c     fjd      - fraction of the julian day.                                    
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
c    *******************************************************************        
c    *                           c o m p j d                           *        
c    *    statement blocked by ray orzol                               *        
c    *******************************************************************        
c                                                                               
                               d i m e n s i o n                                
     1   ndm(12)                                                                
c                                                                               
                                   d a t a                                      
     1   jdor/2415019/,                                                         
     2   jyr19/1900/                                                            
c                                                                               
                                   d a t a                                      
     1   ndm/0,31,59,90,120,151,181,212,243,273,304,334/                        
c    *******************************************************************        
c     computes julian day and fraction from year, month, day and time ut        
c     accurate only between march 1, 1900 and february 28, 2100                 
c     based on julian calendar corrected to correspond to gregorian             
c        calendar during this period                                            
c    jdor=jd of december 30, 1899 at 12 hours ut                                
c    *******************************************************************        
      jd=jdor                                                                   
      jyrm9=jyr-jyr19                                                           
      lp=jyrm9/4                                                                
      if(lp.le.0) go to 4                                                       
      jd=jd+1461*lp                                                             
  4   ny=jyrm9-4*lp                                                             
      ic=0                                                                      
      if(ny.gt.0) go to 5                                                       
      if(jmnth.gt.2) ic=1                                                       
      go to 6                                                                   
  5   jd=jd+365*ny+1                                                            
  6   jd=jd+ndm(jmnth)+jday+ic                                                  
      if(jhr.ge.12) go to 7                                                     
      jd=jd-1                                                                   
      fjd=.5e0+.041666667e0*float(jhr)+.00069444444e0*float(jmn)             
      return                                                                    
  7   fjd=.041666667e0*float(jhr-12)+.00069444444e0*float(jmn)                
      return                                                                    
      end                                                                       
