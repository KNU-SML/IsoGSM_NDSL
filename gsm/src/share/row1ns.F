#define LATG latg
#define LONF lonf
      subroutine row1ns(a)                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    row1ns      pair northern and southern latitudes.              
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-12                     
c                                                                               
c abstract: pair northern and southern latitudes of gaussian grid.              
c                                                                               
c program history log:                                                          
c   88-04-12  joseph sela                                                       
c                                                                               
c usage:    call row1ns (a)                                                     
c   input argument list:                                                        
c     a        - array of gaussian grid with separated                          
c                northern and southern latitudes.                               
c                input array is overwritten by output array.                    
c                                                                               
c   output argument list:                                                       
c     a        - array of gaussian grid with paired                             
c                northern and southern latitudes.                               
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension a(lonf_,latg_)                                                  
      dimension b(lonf_),jdone(latg_)                                           
      do 1 k=1,latg_                                                            
      jdone(k)=0                                                                
 1    continue                                                                  
      jsave=2                                                                   
 3    continue                                                                  
      if(mod(jsave,2).eq.0) then                                                
        jget=latg_+1-jsave/2                                                    
      else                                                                      
        jget=(jsave+1)/2                                                        
      endif                                                                     
      if(jget.eq.jsave) then                                                    
        jdone(jsave)=1                                                          
        go to 35                                                                
      endif                                                                     
      do 5  j=1,lonf_                                                           
      b(j)=a(j,jsave)                                                           
      a(j,jsave)=a(j,jget)                                                      
 5    continue                                                                  
      jdone(jsave)=1                                                            
10    continue                                                                  
      jput=jget                                                                 
      if(mod(jput,2).eq.0) then                                                 
        jget=latg_+1-jput/2                                                     
      else                                                                      
        jget=(jput+1)/2                                                         
      endif                                                                     
      if(jget.eq.jsave) go to 20                                                
      do 15  j=1,lonf_                                                          
      a(j,jput)=a(j,jget)                                                       
15    continue                                                                  
      jdone(jput)=1                                                             
      go to 10                                                                  
20    continue                                                                  
      do 30  j=1,lonf_                                                          
      a(j,jput)=b(j)                                                            
30    continue                                                                  
      jdone(jput)=1                                                             
35    continue                                                                  
      do 40 k=jsave,latg_                                                       
      if(jdone(k).eq.0) then                                                    
        jsave=k                                                                 
        go to 3                                                                 
      endif                                                                     
40    continue                                                                  
      return                                                                    
      end                                                                       
