#define LONF lonf
      subroutine symasa(a,b,n)                                                  
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    symns       forms sum and difference of 2 fouriers.            
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-08                     
c                                                                               
c abstract: forms sum and difference of 2 fouriers.                             
c                                                                               
c program history log:                                                          
c   88-04-08  joseph sela                                                       
c                                                                               
c usage:    call symns (a, b, n)                                                
c   input argument list:                                                        
c     a        - let xxx be array elements a(1,k) ... a(lonf_,k).               
c                let yyy be array elements                                      
c                a(lonf_+1,k) ... a(lonf2_,k).                                  
c                on input xxx is northern hemisphere fouriers.                  
c                on input yyy is southern hemisphere fouriers.                  
c     n        - number of fields of arrays a and b.                            
c                                                                               
c   output argument list:                                                       
c     a        - on output array elements xxx are the sum of                    
c                the northern and southern hemisphere fouriers,                 
c                i.e., the multipliers of symmetric pln.                        
c     b        - on output array b is the difference of                         
c                the northern and southern hemisphere fouriers,                 
c                i.e., the multipliers of antisymmetric pln.                    
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      dimension a(lonf_,2,n), b(lonf_,n)                                        
cc                                                                              
cc    on input a corresp. to n. hemis. fouriers                                 
cc    on input b corresp. to s. hemis. fouriers                                 
cc    on output a corresp. to multipliers of sym. pln                           
cc    on output b corresp. to multipliers of antisym. pln                       
cc                                                                              
      do 10 k=1,n                                                               
      do 200 i=1,lonf_                                                          
      b(i,k)   = a(i,1,k) - a(i,2,k)                                            
  200 continue                                                                  
      do 400 i=1,lonf_                                                          
      a(i,1,k) = a(i,1,k) + a(i,2,k)                                            
  400 continue                                                                  
   10 continue                                                                  
      return                                                                    
      end                                                                       
