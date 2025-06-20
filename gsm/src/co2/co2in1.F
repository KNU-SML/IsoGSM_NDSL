#define LEVS levs
      subroutine co2in1(t20,t21,t66,iq)                                         
c    co2in1=co2ins for method 1                                                 
c     *********************************************************                 
c       save data on permanent data set denoted by co222 ***                    
c          ..... k.campana   march 1988,october 1988                            
c          ..... k.campana   december 88 cleaned up for launcher                
c     *********************************************************                 
#include <paramodel.h>
      parameter (l=levs_,lp1=l+1)                                               
      dimension t20(lp1,lp1,3),t21(lp1,lp1,3),t66(l,6)                          
      dimension dcdt8(lp1,lp1),dcdt10(lp1,lp1),co2po(lp1,lp1),                  
     * co2800(lp1,lp1),co2po1(lp1,lp1),co2801(lp1,lp1),co2po2(lp1,lp1),         
     * co2802(lp1,lp1),n(lp1),d2ct8(lp1,lp1),d2ct10(lp1,lp1)                    
      itin=20                                                                   
      itin1=21                                                                  
co222 latest code has iq=1                                                      
ccc         iq=4                                                                
1011  format (4f20.14)                                                          
ccc        read (itin,1011) ((co2po(i,j),i=1,lp1),j=1,lp1)                      
ccc        read (itin1,1011) ((co2800(i,j),i=1,lp1),j=1,lp1)                    
ccc        read (itin,1011) ((co2po1(i,j),i=1,lp1),j=1,lp1)                     
ccc        read (itin1,1011) ((co2801(i,j),i=1,lp1),j=1,lp1)                    
ccc        read (itin,1011) ((co2po2(i,j),i=1,lp1),j=1,lp1)                     
ccc        read (itin1,1011) ((co2802(i,j),i=1,lp1),j=1,lp1)                    
      do 300 j=1,lp1                                                            
        do 300 i=1,lp1                                                          
          co2po(i,j) = t20(i,j,1)                                               
cnov89                                                                          
          if (iq.eq.5) go to 300                                                
cnov89                                                                          
          co2po1(i,j) = t20(i,j,2)                                              
          co2po2(i,j) = t20(i,j,3)                                              
  300 continue                                                                  
      do 301 j=1,lp1                                                            
        do 301 i=1,lp1                                                          
          co2800(i,j) = t21(i,j,1)                                              
cnov89                                                                          
          if (iq.eq.5) go to 301                                                
cnov89                                                                          
          co2801(i,j) = t21(i,j,2)                                              
          co2802(i,j) = t21(i,j,3)                                              
  301 continue                                                                  
c***the following code is rewritten so that the radiative bands                 
c   are:                                                                        
c        iq=1    560-800     (consol.=490-850)                                  
c        iq=2    560-670     (consol.=490-670)                                  
c        iq=3    670-800     (consol.=670-850)                                  
c        iq=4    560-760 (original code)   (consol.=490-850)                    
cnov89                                                                          
c        iq=5   2270-2380    (consol.=2270-2380)                                
cnov89                                                                          
c  the following loop obtains transmission functions for bands                  
c  used in radiative model calculations,with the equivalent                     
c  widths kept from the original consolidated co2 tf s.                         
      if (iq.eq.1) then                                                         
         c1=1.5                                                                 
         c2=0.5                                                                 
      endif                                                                     
      if (iq.eq.2) then                                                         
        c1=18./11.                                                              
        c2=7./11.                                                               
      endif                                                                     
      if (iq.eq.3) then                                                         
        c1=18./13.                                                              
        c2=5./13.                                                               
      endif                                                                     
      if (iq.eq.4) then                                                         
        c1=1.8                                                                  
        c2=0.8                                                                  
      endif                                                                     
cnov89                                                                          
      if (iq.eq.5) then                                                         
        c1=1.0                                                                  
        c2=0.0                                                                  
      endif                                                                     
cnov89                                                                          
      do 1021 i=1,lp1                                                           
      do 1021 j=1,lp1                                                           
      co2po(j,i)=c1*co2po(j,i)-c2                                               
      co2800(j,i)=c1*co2800(j,i)-c2                                             
cnov89                                                                          
      if (iq.eq.5) go to 1021                                                   
cnov89                                                                          
      co2po1(j,i)=c1*co2po1(j,i)-c2                                             
      co2801(j,i)=c1*co2801(j,i)-c2                                             
      co2po2(j,i)=c1*co2po2(j,i)-c2                                             
      co2802(j,i)=c1*co2802(j,i)-c2                                             
1021  continue                                                                  
cnov89                                                                          
      if (iq.ge.1.and.iq.le.4) then                                             
cnov89                                                                          
      do 1 j=1,lp1                                                              
      do 1 i=1,lp1                                                              
      dcdt8(i,j)=.02*(co2801(i,j)-co2802(i,j))*100.                             
      dcdt10(i,j)=.02*(co2po1(i,j)-co2po2(i,j))*100.                            
      d2ct8(i,j)=.0016*(co2801(i,j)+co2802(i,j)-2.*co2800(i,j))*1000.           
      d2ct10(i,j)=.0016*(co2po1(i,j)+co2po2(i,j)-2.*co2po(i,j))*1000.           
1     continue                                                                  
cnov89                                                                          
      endif                                                                     
cnov89                                                                          
co222 *********************************************************                 
ccc          rewind 66                                                          
c        save cdtm51,co2m51,c2dm51,cdtm58,co2m58,c2dm58..on tempo file          
ccc          write (66) (dcdt10(i,i+1),i=1,l)                                   
ccc          write (66) (co2po(i,i+1),i=1,l)                                    
ccc          write (66) (d2ct10(i,i+1),i=1,l)                                   
ccc          write (66) (dcdt8(i,i+1),i=1,l)                                    
ccc          write (66) (co2800(i,i+1),i=1,l)                                   
ccc          write (66) (d2ct8(i,i+1),i=1,l)                                    
ccc          rewind 66                                                          
co222 *********************************************************                 
      do 400 i=1,l                                                              
        t66(i,2) = co2po(i,i+1)                                                 
        t66(i,5) = co2800(i,i+1)                                                
cnov89                                                                          
        if (iq.eq.5) go to 400                                                  
cnov89                                                                          
        t66(i,1) = dcdt10(i,i+1)                                                
        t66(i,3) = d2ct10(i,i+1)                                                
        t66(i,4) = dcdt8(i,i+1)                                                 
        t66(i,6) = d2ct8(i,i+1)                                                 
  400 continue                                                                  
      return                                                                    
      end                                                                       
