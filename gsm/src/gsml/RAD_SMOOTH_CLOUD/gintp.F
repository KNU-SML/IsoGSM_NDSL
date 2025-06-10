      subroutine gintp(iin,jtwidl,jin,iout,                                     
     1                 inslat,wgtlat,                               
     2                 ctt,cldt,tt,sum,nn,ltwidl,latrd1)                        
      dimension ctt(iin,jtwidl)                                                 
      dimension cldt(iout)                                                      
      dimension tt(iout,4),sum(iout)                                            
      dimension nn(iout)                                                        
c....    for top temp just take average of surrounding                          
c..       non-zero points (these are the cloud-filled ones)....                 
c..       nn will be number of surrounding pts with cld (gt zero)               
      lonf=iin/2
      if (inslat.lt.0) go to 600                                                
      inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1              
      inth1 = mod(inth,jtwidl) + 1                                              
      if (inslat.eq.jin) go to 105                                              
      do 100 i=1,iout                                                           
        ileft=i
        irght=i+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        tt(i,1) = ctt(ileft,inth)                                            
        tt(i,2) = ctt(ileft,inth1)                                           
        tt(i,3) = ctt(irght,inth)                                            
        tt(i,4) = ctt(irght,inth1)                                           
  100 continue                                                                  
      go to 130                                                                 
  105 do 110 i=1,iout                                                           
        ileft=i
        irght=i+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        tt(i,1) = ctt(ileft,inth)                                            
        tt(i,3) = ctt(irght,inth)                                            
  110 continue                                                                  
      iout2 = iout / 2                                                          
      do 120 i=1,iout2                                                          
        ileft=i
        irght=ileft+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        ileft2=i+iout2
        irght2=ileft2+1
        if(mod(ileft2,lonf).eq.0) irght2=irght2-lonf
        tt(i      ,2) = ctt(ileft2,inth)                                      
        tt(i+iout2,2) = ctt(ileft ,inth)                                      
        tt(i      ,4) = ctt(irght2,inth)                                      
        tt(i+iout2,4) = ctt(irght ,inth)                                      
  120 continue                                                                  
c---      nn will be number of surrounding pts with cld (gt zero)               
  130 do 10 i=1,iout                                                            
        nn(i) = 0                                                               
   10 continue                                                                  
      do 12 i=1,iout                                                            
        sum(i) = 0.e0                                                          
   12 continue                                                                  
      do 150 kpt=1,4                                                            
       do 14 i=1,iout                                                           
        if (tt(i,kpt).gt.0.e0) then                                            
          nn(i) = nn(i) + 1                                                     
          sum(i) = sum(i) + tt(i,kpt)                                           
        endif                                                                   
   14  continue                                                                 
  150  continue                                                                 
      do 16 i=1,iout                                                            
        if (nn(i).lt.1) then                                                    
          cldt(i) = 0.e0                                                       
        else                                                                    
          cldt(i) = sum(i) / nn(i)                                              
        end if                                                                  
   16 continue                                                                  
      return                                                                    
c--- polar region                                                               
  600 continue                                                                  
      ja = iabs(inslat)                                                         
      do 200 i=1,iout                                                           
        ileft=i
        irght=ileft+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        tt(i,1) = ctt(ileft,ja)                                              
        tt(i,2) = ctt(irght,ja)                                              
  200 continue                                                                  
c---      nn will be number of surrounding pts with cld (gt zero)               
      do 20 i=1,iout                                                            
        nn(i) = 0                                                               
   20 continue                                                                  
      do 22 i=1,iout                                                            
       sum(i) = 0.e0                                                           
   22 continue                                                                  
      do 202 kpt=1,2                                                            
        do 24 i=1,iout                                                          
          if (tt(i,kpt).gt.0.e0) then                                          
            nn(i) = nn(i) + 1                                                   
            sum(i) = sum(i) + tt(i,kpt)                                         
          endif                                                                 
   24   continue                                                                
  202 continue                                                                  
      do 26 i=1,iout                                                            
        if (nn(i).lt.1) then                                                    
          cldt(i) = 0.e0                                                       
        else                                                                    
          cldt(i) = sum(i) / nn(i)                                              
        end if                                                                  
   26 continue                                                                  
      return                                                                    
      end                                                                       
