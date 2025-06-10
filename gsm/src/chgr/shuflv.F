#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine shuflv(tri, rho, iswtch)                                       
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    shuflv      vector triangle <--> rhomboid shuffle.             
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-21                     
c                                                                               
c abstract: vector triangle <--> rhomboid shuffle.                              
c                                                                               
c program history log:                                                          
c   88-04-21  joseph sela                                                       
c                                                                               
c usage:    call shuflv (tri, rho, iswitch)                                     
c   input argument list:                                                        
c     tri      - triangle.                                                      
c                when iswitch=-1, triangle tri is stored in                     
c                the lower triangle of rhomboid rho, and                        
c                the upper triangle of rhomboid rho is set to zeros.            
c     rho      - rhomboid.                                                      
c                when iswitch=1, the lower triangle of rhomboid rho             
c                is stored in triangle tri.                                     
c     iswitch  - integer switch.                                                
c                shufle must be called first with iswitch=0                     
c                to initialize bit vector bitv.                                 
c                                                                               
c   output argument list:                                                       
c     tri      - triangle.                                                      
c                when iswitch=1, the lower triangle of rhomboid rho             
c                is stored in triangle tri.                                     
c     rho      - rhomboid.                                                      
c                when iswitch=-1, triangle tri is stored in                     
c                the lower triangle of rhomboid rho, and                        
c                the upper triangle of rhomboid rho is set to zeros.            
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
      save                                                                      
#include <paramter.h>
c                                                                               
      parameter(idimt=idim*2,jdimhf=jdim/2,ijdim=idim*jdim,                     
     1          mdim=(mwave+1)*(mwave+2),                                       
     3          mdimv=(mwave+1)*(mwave+4),                                      
     4          mwvp2=(mwave+1)*2,mwave2=mwave+2,mwavep=mwave+1,                
     5          mdimhf=mdim/2,mdimvh=mdimv/2,                                   
     6          kdimp=kdim+1,kdimm=kdim-1)                                      
      dimension tri( mdimv ), rho( mwvp2 , mwave2 )                             
c                                                                               
c     when iswtch=1, the lower triangle of rhomboid rho                         
c     is stored in triangle tri.                                                
c                                                                               
c     when iswtch=-1, triangle tri is stored in                                 
c     the lower triangle of rhomboid rho, and                                   
c     the upper triangle of rhomboid rho is set to zeros.                       
c                                                                               
      if (iswtch) 200, 400, 600                                                 
c                                                                               
  200 continue                                                                  
      do 340 k=1, mwave2                                                        
      do 320 i=1, mwvp2                                                         
             rho(i,k) = 0.0                                                     
  320 continue                                                                  
  340 continue                                                                  
      do 350 i=1, mwvp2                                                         
             rho(i,1) = tri(i)                                                  
  350 continue                                                                  
      iplus =  mwvp2                                                            
      len   =  mwvp2                                                            
      do 380 k=2, mwave2                                                        
      do 360 i=1,len                                                            
             rho(i,k) = tri(i+iplus)                                            
  360 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  380 continue                                                                  
      return                                                                    
c                                                                               
  400 continue                                                                  
      return                                                                    
c                                                                               
  600 continue                                                                  
      do 750 i=1, mwvp2                                                         
             tri(i) = rho(i,1)                                                  
  750 continue                                                                  
      iplus =  mwvp2                                                            
      len   =  mwvp2                                                            
      do 780 k=2, mwave2                                                        
      do 760 i=1,len                                                            
             tri(i+iplus) = rho(i,k)                                            
  760 continue                                                                  
      iplus = iplus + len                                                       
      len = len - 2                                                             
  780 continue                                                                  
      return                                                                    
      end                                                                       
