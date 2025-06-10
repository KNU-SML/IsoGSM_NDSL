#define IGRD igrd
#define IGRD1 igrd1
#define JGRD jgrd
#define JGRD1 jgrd1
#define JGRD12 jgrd12
#define G2G
      subroutine shalfo(a,km)                                                   
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:  shalfo                                                           
c   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06          
c                                                                               
c abstract:  transform model grid (i,j) to regular grid (2i,j/2) with           
c            north and south grid combined together as global model for         
c            regional model computation.                                        
c            the model grid is                                                  
c               (1,1),(2,1),.......(i,1),(1,j  ),(2,j  ),...,(i,j  )            
c               (1,2),(2,2),......,(i,2),(1,j-1),(2,j-1),...,(i,j-1)            
c                ...................................................            
c               (1,j/2),(2,j/2),..,(i,j/2),(1,j-j/2+1),...,(i,j-j/2+1)          
c            the regular grid is                                                
c               (1,1),(2,1),.......(i,1)                                        
c               (1,2),(2,2),.......(i,2)                                        
c               ........................                                        
c               (1,j),(2,j),.......(i,j)                                        
c                                                                               
c program history log:                                                          
c                                                                               
c usage:    call shalfo(a,k)                                                    
c   input argument list:                                                        
c     a 	- model grid with dimension (i,j,k)                                    
c     k 	- the third dimension of a                                             
c                                                                               
c   output argument list:                                                       
c     a 	- regular output grid with dimension (i,j,k)                           
c                                                                               
c   input files: none                                                           
c                                                                               
c   output files: none                                                          
c                                                                               
c   subprograms called: none                                                    
c                                                                               
c   remark: none                                                                
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
#include <paramodel.h>
#ifndef G2G                                                                                
c                                                                               
c arrage j row for output (regional regular grid)                               
c                                                                               
      dimension a(igrd1_,jgrd1_,km)                                             
      dimension tmp(igrd1_,jgrd1_)                                              
c                                                                               
      do 10 k=1,km                                                              
      do 20 j=1,jgrd12_                                                         
      jj = 2*j-1                                                                
      js = j                                                                    
      jn = jgrd1_+1 - j                                                         
      do 30 i=1,igrd1_                                                          
      tmp(i,js) = a(i,jj  ,k)                                                   
      tmp(i,jn) = a(i,jj+1,k)                                                   
 30   continue                                                                  
 20   continue                                                                  
      do 40 j=1,jgrd1_                                                          
      do 40 i=1,igrd1_                                                          
      a(i,j,k) = tmp(i,j)                                                       
 40   continue                                                                  
 10   continue                                                                  
c                                                                               
#endif
      return                                                                    
      end                                                                       
