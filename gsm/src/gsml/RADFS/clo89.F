#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp

      subroutine clo89(ipts,cldfac,camt,nclds,kbtm,ktop)                             
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c                                                                               
#include <paramodel.h>
#include <hcon.h>
c     subroutine clo88 computes cloud transmission functions for the            
c  longwave code,using code written by bert katz (301-763-8161).                
c  and modified by dan schwarzkopf in december,1988.                            
c                inputs:          (common block)                                
c      camt,ktop,kbtm,nclds         radisw                                      
c                output:                                                        
c      cldfac                       cldcom                                      
c                                                                               
c          called by:      radmn or model routine                               
c          calls    :                                                           
c                                                                               
#include <rdparm.h>
c                                                                               
      dimension nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)                       
      dimension camt(imbx,lp1),cldfac(imbx,lp1,lp1)                             
      dimension cldrow(lp1)                                                     
c                                                                               
c      dimension cldfip(lp1,lp1)                                                
       dimension cldipt(lp1,lp1,NVECT)                                         
c                                                                               
      do 1 iq=1,ipts
      itop=iq
      jtop=1                                                            
      do 11 ip=1,jtop                                                           
      ir=iq+ip-1                                                                
      if (nclds(ir).eq.0) then                                                  
        do 25 j=1,lp1                                                           
        do 25 i=1,lp1                                                           
        cldipt(i,j,ip)=1.                                                       
25      continue                                                                
      endif                                                                     
      if (nclds(ir).ge.1) then                                                  
          xcld=1.-camt(ir,2)                                                    
           k1=ktop(ir,2)+1                                                      
           k2=kbtm(ir,2)                                                        
          do 27 j=1,lp1                                                         
              cldrow(j)=1.                                                      
27        continue                                                              
          do 29 j=1,k2                                                          
              cldrow(j)=xcld                                                    
29        continue                                                              
          kb=max(k1,k2+1)                                                       
          do 33 k=kb,lp1                                                        
          do 33 kp=1,lp1                                                        
               cldipt(kp,k,ip)=cldrow(kp)                                       
33        continue                                                              
          do 37 j=1,lp1                                                         
              cldrow(j)=1.                                                      
37        continue                                                              
          do 39 j=k1,lp1                                                        
              cldrow(j)=xcld                                                    
39        continue                                                              
          kt=min(k1-1,k2)                                                       
          do 43 k=1,kt                                                          
          do 43 kp=1,lp1                                                        
              cldipt(kp,k,ip)=cldrow(kp)                                        
43        continue                                                              
          if(k2+1.le.k1-1) then                                                 
            do 31 j=k2+1,k1-1                                                   
            do 31 i=1,lp1                                                       
                cldipt(i,j,ip)=1.                                               
31          continue                                                            
          else if(k1.le.k2) then                                                
            do 32 j=k1,k2                                                       
            do 32 i=1,lp1                                                       
                cldipt(i,j,ip)=xcld                                             
32          continue                                                            
          endif                                                                 
      endif                                                                     
      if (nclds(ir).ge.2) then                                                  
        do 21 nc=2,nclds(ir)                                                    
          xcld=1.-camt(ir,nc+1)                                                 
           k1=ktop(ir,nc+1)+1                                                   
           k2=kbtm(ir,nc+1)                                                     
          do 47 j=1,lp1                                                         
              cldrow(j)=1.                                                      
47        continue                                                              
          do 49 j=1,k2                                                          
              cldrow(j)=xcld                                                    
49        continue                                                              
          kb=max(k1,k2+1)                                                       
          do 53 k=kb,lp1                                                        
          do 53 kp=1,lp1                                                        
               cldipt(kp,k,ip)=cldipt(kp,k,ip)*cldrow(kp)                       
c              cldfip(kp,k)=cldrow(kp)                                          
53        continue                                                              
          do 57 j=1,lp1                                                         
              cldrow(j)=1.                                                      
57        continue                                                              
          do 59 j=k1,lp1                                                        
              cldrow(j)=xcld                                                    
59        continue                                                              
          kt=min(k1-1,k2)                                                       
          do 63 k=1,kt                                                          
          do 63 kp=1,lp1                                                        
              cldipt(kp,k,ip)=cldipt(kp,k,ip)*cldrow(kp)                        
c             cldfip(kp,k)=cldrow(kp)                                           
63        continue                                                              
c         if(k2+1.le.k1-1) then                                                 
c           do 51 j=k2+1,k1-1                                                   
c           do 51 i=1,lp1                                                       
c               cldipt(i,j,ip)=1.                                               
c51         continue                                                            
          if(k1.le.k2) then                                                     
            do 52 j=k1,k2                                                       
            do 52 i=1,lp1                                                       
                cldipt(i,j,ip)=cldipt(i,j,ip)*xcld                              
52          continue                                                            
          endif                                                                 
c           do 65 j=1,lp1                                                       
c           do 65 i=1,lp1                                                       
c         cldipt(i,j,ip)=cldipt(i,j,ip)*cldfip(i,j)                             
c65       continue                                                              
21        continue                                                              
      endif                                                                     
11    continue                                                                  
      do 71 j=1,lp1                                                             
      do 71 i=1,lp1                                                             
      do 71 ip=1,jtop                                                           
      ir=iq+ip-1                                                                
      cldfac(ir,i,j)=cldipt(i,j,ip)                                             
71    continue                                                                  
1     continue                                                                  
      return                                                                    
      end                                                                       
