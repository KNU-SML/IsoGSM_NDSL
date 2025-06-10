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
      subroutine cmpind                                                         
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    cmpind      compute indices for spectral transpose.            
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: computes indices to transpose spectral coefficients                 
c           between ibm order (data run up the columns)                         
c           and cyber order (data run up the diagonals).                        
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell  docblock written (prehistorical program)            
c                                                                               
c usage:    call cmpind                                                         
c   commons:                                                                    
c     comind   - indices indxnn for transposing ibm to cyber order,             
c                indices indxmm for transposing cyber to ibm order.             
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  
c                                                                               
c$$$                                                                            
c                                                                               
#include <paramter.h>
c                                                                               
      parameter(mwavep=mwave+1,                                                 
     1          mdim=(mwave+1)*(mwave+2))                                       
c                                                                               
      common/comind/ indxnn(mdim),indxmm(mdim)                                  
c                                                                               
c  indxnn(mdim) :  1-d index of converting input form spher coeff array         
c                  to transposed form array                                     
c  indxmm(mdim) :  1-d index of converting transposed form spher coeff          
c                  array to input form spherical coeff array                    
c                                                                               
      l=0                                                                       
      do 10 m=1,mwavep                                                          
      nend=mwavep-m+1                                                           
      do 10 nn=1,nend                                                           
      n=nn+m-1                                                                  
      l=l+2                                                                     
      indx=(mwavep*(n-m)-(n-m)*(n-m-1)/2+m)*2-1                                 
      indxnn(l-1)=indx                                                          
      indxnn(l  )=indx+1                                                        
   10 continue                                                                  
c                                                                               
      l=0                                                                       
      do 20 nn=1,mwavep                                                         
      lln=mwavep-nn+1                                                           
      do 20 ll=1,lln                                                            
      n=ll+nn-1                                                                 
      m=ll                                                                      
      indx=(m*mwavep-(mwavep-n)-(m-1)*m/2)*2-1                                  
      l=l+2                                                                     
      indxmm(l-1)=indx                                                          
      indxmm(l  )=indx+1                                                        
   20 continue                                                                  
      return                                                                    
c                                                                               
      end                                                                       
