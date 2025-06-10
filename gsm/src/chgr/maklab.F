      subroutine maklab(lab,idate,fhour,ifin)                                   
      save                                                                      
c$$$  subprogram documentation block                                            
#include <paramodel.h>
c                .      .    .                                       .          
c subprogram:    maklab      sets beginning of office note 85 label.            
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-29                     
c                                                                               
c abstract: sets first 8 bytes of office note 85 label.                         
c                                                                               
c program history log:                                                          
c   88-04-29  joseph sela                                                       
c                                                                               
c usage:    call maklab (lab, idate, fhour, ifin)                               
c   input argument list:                                                        
c     lab      - office note 85 label.                                          
c     idate    - idate(1)=initial hour (gmt) of forecast.                       
c                idate(2)=month (1-12).                                         
c                idate(3)=day of the month.                                     
c                idate(4)=year.                                                 
c     fhour    - forecast hour.                                                 
c     ifin     - integer switch.                                                
c                when ifin.gt.0,                                                
c                first 8 bytes of lab are set to                                
c                'siginib2', 'sigin6b2', 'sigge6b2', 'siggesb2',                
c                depending on values of fhour and idate(1).                     
c                when ifin.le.0,                                                
c                first 8 bytes of lab are set to                                
c                'smsxxxb2', where xxx is the forecast hour, fhour.             
c                                                                               
c   output argument list:                                                       
c     lab      - office note 85 label with first 8 bytes set.                   
c                                                                               
c   output files:                                                               
c     output   - print file.                                                    
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
c                                                                               
      dimension idigts(3)                                                       
      dimension idate (4)                                                       
c                                                                               
      character*4 lab (8)                                                       
      character*4 lchars                                                        
c                                                                               
      ihour=fhour+0.5                                                           
      write(lchars,103) ihour                                                   
 103  format (i4)                                                               
      read(lchars,105) idigts                                                   
 105  format (1x, 3i1)                                                          
      write(lab(1),107) idigts(1)
 107  format (3x,i1)                                                      
      write(lab(2),108) idigts(2),idigts(3)
 108  format (2i1,'b2')                                                      
c                                                                               
      return                                                                    
      end                                                                       
