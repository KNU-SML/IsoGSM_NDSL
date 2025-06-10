      subroutine rctrns(itap1,ratstd,ratsm,ratio)                               
c      ratstd=value of higher std co2 concentration                             
c      ratsm=value of lower std co2 concentration                               
c      ratio=actual co2 concentration                                           
c      the 3 above quantities are in units of 330 ppmv.                         
      common/input/p1,p2,trnslo,ia,ja,n                                         
      common/press/pa(109)                                                      
      common/tran/transa(109,109)                                               
      dimension trns1(109,109),trns2(109,109)                                   
100   format (4f20.14)                                                          
c   read in tfs of lower std co2 concentration                                  
      read (itap1,100) ((trns1(i,j),i=1,109),j=1,109)                           
      itap2=itap1+1                                                             
c   read in tfs of higher std co2 concentration                                 
      read (itap2,100) ((transa(i,j),i=1,109),j=1,109)                          
      call coeint(ratstd)                                                       
      do 401 i=1,109                                                            
      do 401 j=1,i                                                              
      if (j.eq.i) go to 401                                                     
c  using higher co2 concentration,compute 1st guess co2 tfs for                 
c  actual co2 concentration.                                                    
      p2=(ratio+ratstd)*pa(i)/(2.*ratstd)  +                                    
     1   (ratstd-ratio)*pa(j)/(2.*ratstd)                                       
      p1=(ratstd-ratio)*pa(i)/(2.*ratstd)  +                                    
     1   (ratio+ratstd)*pa(j)/(2.*ratstd)                                       
      call sintr2                                                               
      trnspr=trnslo                                                             
c  using higher co2 concentration,compute 1st guess co2 tfs for                 
c  lower std co2 concentration                                                  
      p2=(ratsm+ratstd)*pa(i)/(2.*ratstd)  +                                    
     1   (ratstd-ratsm)*pa(j)/(2.*ratstd)                                       
      p1=(ratstd-ratsm)*pa(i)/(2.*ratstd)  +                                    
     1   (ratsm+ratstd)*pa(j)/(2.*ratstd)                                       
      call sintr2                                                               
      trnspm=trnslo                                                             
c  compute tfs for co2 concentration given by (ratio).                          
c   store temporarily in (trns2)                                                
      trns2(j,i)=trnspr+(ratstd-ratio)*(trns1(j,i)-                             
     1 trnspm)/(ratstd-ratsm)                                                   
      trns2(i,j)=trns2(j,i)                                                     
c we now can overwrite (trns1) and store in (trns1) the 1st guess               
c  co2 tfs for lower std co2 concentration                                      
      trns1(j,i)=trnslo                                                         
      trns1(i,j)=trnslo                                                         
401   continue                                                                  
c  set diagonal values of co2 tfs to unity                                      
      do 402 i=1,109                                                            
      trns1(i,i)=1.0                                                            
      trns2(i,i)=1.0                                                            
402   continue                                                                  
c  now output the computed co2 tfs for (ratio) co2 conc. in (transa)            
      do 405 i=1,109                                                            
      do 405 j=1,109                                                            
      transa(j,i)=trns2(j,i)                                                    
405   continue                                                                  
      return                                                                    
      end                                                                       
