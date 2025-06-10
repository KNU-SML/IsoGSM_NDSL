#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF lonf
#define NPES npes
      subroutine rowsep(a)                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    rowsep      separate northern and southern latitudes.          
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-12                     
c                                                                               
c abstract: separate northern and southern latitudes of gaussian grid.          
c                                                                               
c program history log:                                                          
c   88-04-12  joseph sela                                                       
c                                                                               
c usage:    call rowsep (a)                                                     
c   input argument list:                                                        
c     a        - array of gaussian grid with paired                             
c                northern and southern latitudes.                               
c                input array is overwritten by output array.                    
c                                                                               
c   output argument list:                                                       
c     a        - array of gaussian grid with separated                          
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
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c
#ifdef REDUCE_GRID
      lonf2=lonf_*2
      latg2=latg_/2
      call reduceintp(a,lonfd,lonf2,latg2)
#endif
c
      do 1 k=1,latg_                                                            
      jdone(k)=0                                                                
 1    continue                                                                  
      jsave=2                                                                   
 3    continue                                                                  
      if(jsave.le.latg2_) then                                                  
        jget=2*jsave-1                                                          
      else                                                                      
        jget=2*(latg_+1-jsave)                                                  
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
      if(jput.le.latg2_) then                                                   
        jget=2*jput-1                                                           
      else                                                                      
        jget=2*(latg_+1-jput)                                                   
      endif                                                                     
      if(jget.eq.jsave) go to 20                                                
      do 15 j=1,lonf_                                                           
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
c
      return                                                                    
      end                                                                       
