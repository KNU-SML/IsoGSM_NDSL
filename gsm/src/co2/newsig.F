#define CP cp
#define JCAP jcap
#define LATG latg
#define LEVS levs
#define LONF lonf
#define LSOIL lsoil
#define RD rd
      subroutine newsig(ci, si, del, sl, cl, rpi)                               
c                                                                               
#include <paramodel.h>
       parameter(                                                               
     *           latdyn=latg_,                                                  
     *           londyn=lonf_,                                                  
     *           lonf=londyn,                                                   
     *           latg=latdyn)                                                   
c                                                                               
      save                                                                      
c ------------------********************-----------------------                 
c this routine sets coordinates for levels and use phillips method              
c to get layers. ( can be stretched levels or not )                             
c  --------------------------------------------------------------               
c                                                                               
#include <paramco2.h>
c                                                                               
      parameter(kdimp=kdim+1,kdimm=kdim-1)                                      
c                                                                               
      dimension ci(1),si(1),del(1),sl(1),cl(1),rpi(1)                           
c                                                                               
      dimension delsig(kdim)                                                    
c                                                                               
      real rk,rk1,rkinv                                                         
c                                                                               
c  model sigma data in array delmdl of length levmdl                            
#include <modlsigs.h>
c                                                                               
      data ifp/0/                                                               
c***********************************************************************        
c  in case of shocking parameter discord snafu---                               
      print *,'levmdl,kdim=',levmdl,kdim                                        
      if (levmdl.ne.kdim)  stop 18                                              
c                                                                               
      do 1 k=1,kdim                                                             
    1 del(k)=delmdl(k)                                                          
c                                                                               
      ci(1) = 0.                                                                
      do 2 k=1,kdim                                                             
    2 ci(k+1)=ci(k)+del(k)                                                      
      ci(kdimp)=1.                                                              
c                                                                               
      rk  = rd_/cp_                                                             
      rk1 = rk + 1.                                                             
      rkinv=1./rk                                                               
c                                                                               
      levs=kdim                                                                 
c                                                                               
      do 3 li=1,kdimp                                                           
    3 si(li) = 1. - ci(li)                                                      
c                                                                               
      do 4 le=1,kdim                                                            
      dif = si(le)**rk1 - si(le+1)**rk1                                         
      dif = dif /(rk1*(si(le)-si(le+1)))                                        
      sl(le) = dif**rkinv                                                       
      cl(le) = 1. - sl(le)                                                      
    4 continue                                                                  
c                                                                               
c     compute pi ratios for temp. matrix.                                       
c                                                                               
      do 5 le=1,kdimm                                                           
      rpi(le) = (sl(le+1)/sl(le))                                               
    5 continue                                                                  
      do 15 le=1,kdimm                                                          
      rpi(le) = rpi(le)**rk                                                     
   15 continue                                                                  
c                                                                               
      do 6 le=1,kdimp                                                           
      print 100, le, ci(le), si(le)                                             
100   format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
    6 continue                                                                  
c                                                                               
      print 200                                                                 
  200 format (1h0)                                                              
      do 7 le=1,kdim                                                            
      print 300, le, cl(le), sl(le), del(le)                                    
  300 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,          
     1 'del=', f6.3)                                                            
    7 continue                                                                  
c                                                                               
      print 400, (rpi(le), le=1,kdimm)                                          
  400 format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
c                                                                               
      return                                                                    
      end                                                                       
