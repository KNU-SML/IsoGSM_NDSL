      subroutine cdate(jd,fjd,munth,im,id,iyear,ihr,xmin)                       
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    cdate       computes day,month,yr from julian day              
c   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07                     
c                                                                               
c abstract: this code written at gfdl ....                                      
c   computes month,day,year from julian day.                                    
c   accurate only between march 1, 1900 and february 28, 2100....               
c   based on julian calender corrected to correspond to gregorian               
c   calender during this period.                                                
c                                                                               
c program history log:                                                          
c   77-06-07  robert white,gfdl                                                 
c                                                                               
c usage:    call cdate(jd,fjd,munth,im,id,iyear,ihr,xmin)                       
c   input argument list:                                                        
c     jd       - julian day for current fcst hour.                              
c     fjd      - fraction of the julian day.                                    
c   output argument list:                                                       
c     munth    - month (character).                                             
c     im       - month (integer).                                               
c     id       - day of the month.                                              
c     iyear    - year.                                                          
c     ihr      - hour of the day.                                               
c     xmin     - minute of the hour.                                            
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
c                                                                               
c    *******************************************************************        
c    *                                                                 *        
c    *                            c d a t e                            *        
c    *                                                                 *        
c    *******************************************************************        
c                                                                               
c         statements blocked by robert k. white.......7 june 1977               
c                                                                               
c..... cdate computes month, day, and year from julian day.                     
c.....accurate only between march 1, 1900 and february 28, 2100                 
c.....based on julian calender corrected to correspond to gregorian             
c.....calender during this period                                               
c                                                                               
                              d i m e n s i o n                                 
     1   dy(13),             month(12)                                          
c                                                                               
c                                                                               
                                   d a t a                                      
     1   dy                  /                                                  
     2   0.,                 31.,                59.,                           
     3   90.,                120.,               151.,                          
     4   181.,               212.,               243.,                          
     5   273.,               304.,               334.,                          
     6   365.                /                                                  
c                                                                               
                                   d a t a                                      
     1   month               /                                                  
     2   4hjan.,             4hfeb.,             4hmar.,                        
     3   4hapr.,             4hmay ,             4hjune,                        
     4   4hjuly,             4haug.,             4hsep.,                        
     5   4hoct.,             4hnov.,             4hdec.                         
     6   /                                                                      
c                                                                               
c.....jdor = jd of december 30, 1899 at 12 hours ut                             
c                                                                               
                                   d a t a                                      
     1   jdor                /         2415019             /,                   
     2   iyr                 /         1900                /                    
c                                                                               
c    *******************************************************************        
c                                                                               
      iyear=iyr                                                                 
      nday=jd-jdor                                                              
      if(fjd.ge..5e0) nday=nday+1                                              
 61   if(nday.lt.1462) go to 62                                                 
      nday=nday-1461                                                            
      iyear=iyear+4                                                             
      go to 61                                                                  
 62   ndiy=365                                                                  
      if(mod(iyear,4).eq.0) ndiy=366                                            
      if(nday.le.ndiy) go to 65                                                 
      iyear=iyear+1                                                             
      nday=nday-ndiy                                                            
      go to 62                                                                  
 65   if(nday.gt.int(dy(2))) go to 66                                           
      im=1                                                                      
      id=nday                                                                   
      go to 67                                                                  
 66   if(nday.ne.60) go to 68                                                   
      if(ndiy.eq.365) go to 68                                                  
      im=2                                                                      
      id=29                                                                     
      go to 67                                                                  
 68   if(nday.gt.(int(dy(3))+ndiy-365)) go to 69                                
      im=2                                                                      
      id=nday-31                                                                
      go to 67                                                                  
 69   do 70 i=3,12                                                              
      if(nday.gt.(int(dy(i+1))+ndiy-365)) go to 70                              
      im=i                                                                      
      id=nday-int(dy(i))-ndiy+365                                               
      go to 67                                                                  
 70   continue                                                                  
 67   munth=month(im)                                                           
      hr=24.e0*fjd                                                             
      ihr=hr                                                                    
      xmin=60.e0*(hr-float(ihr))                                               
      ihr=ihr+12                                                                
      if(ihr.ge.24) ihr=ihr-24                                                  
      return                                                                    
      end                                                                       
