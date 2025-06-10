#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LSOILIN lsoilin
#define LTSTP ltstp
#define LVAR lvar
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine transp(a,kmax)                                                 
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    transp      transposes spectral coefficients.                  
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-20                     
c                                                                               
c abstract: transposes rhomboidal spectral coefficients.                        
c                                                                               
c program history log:                                                          
c   88-04-20  joseph sela                                                       
c                                                                               
c usage:    call transp (a, n)                                                  
c   input argument list:                                                        
c     a        - rhomboidal fields of spectral coefficients.                    
c     n        - number of  fields in array a.                                  
c                                                                               
c   output argument list:                                                       
c     a        - transposed fields of spectral coefficients.                    
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
c                                                                               
#include <paramter.h>
c                                                                               
      parameter(idimt=idim*2,jdimhf=jdim/2,                                     
     1          mwavep=mwave+1,mwave2=mwave+2,                                  
     2          mdim=(mwave+1)*(mwave+2),                                       
     3          mdimv=(mwave+1)*(mwave+4),                                      
     4          mwvp2=(mwave+1)*2,                                              
     5          mdimhf=mdim/2,mdimvh=mdimv/2,                                   
     6          kdimp=kdim+1,kdimm=kdim-1)                                      
c                                                                               
      common/comind/ indxnn(mdim),indxmm(mdim)                                  
c                                                                               
      dimension a(mdim,kmax)                                                    
c                                                                               
#include <padchgr.h>
      parameter(npad=nscrch-ntrns)                                              
      common /scratch/ b(mdim),pad(npad)                                        
c                                                                               
      do 10 k=1,kmax                                                            
      do 1 i = 1, mdim                                                          
      b(indxnn(i)) = a(i,k)                                                     
    1 continue                                                                  
      do 2 i = 1, mdim                                                          
      a(i,k)=b(i)                                                               
    2 continue                                                                  
 10   continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
