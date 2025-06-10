#define CP cp
#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
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
#define RD rd
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine setsig(ci, si, del, sl, cl, rpi)                               
#include <paramodel.h>
      save                                                                      
c                                                                               
c modified by hann-ming henry juang for gread  (read sigma file)                
c sigma file contains si and sl, so                                             
c input: si sl                                                                  
c output: ci cl del rpi                                                         
c                                                                               
#include <paramter.h>
#include <parmchgr.h>
c                                                                               
      parameter(kdimpi=kdimi+1,kdimmi=kdimi-1)                                  
c                                                                               
      dimension ci(kdimpi),si(kdimpi),del(kdimi),sl(kdimi),cl(kdimi),           
     &          rpi(kdimmi)                                                     
      real rk,rk1,rkinv                                                         
c                                                                               
      rk  = rd_/cp_                                                             
      rk1 = rk + 1.                                                             
      rkinv=1./rk                                                               
c                                                                               
      do 10 k=1,kdimpi                                                          
      ci(k)=1. - si(k)                                                          
 10   continue                                                                  
c                                                                               
      do 15 k=1,kdimi                                                           
      del(k)=ci(k+1)-ci(k)                                                      
 15   continue                                                                  
c                                                                               
      do 20 k=1,kdimi                                                           
      cl(k)=1. - sl(k)                                                          
 20   continue                                                                  
c                                                                               
c     compute pi ratios for temp. matrix.                                       
c                                                                               
      do 25 le=1,kdimmi                                                         
      rpi(le) = (sl(le+1)/sl(le))                                               
 25   continue                                                                  
      do 30 le=1,kdimmi                                                         
      rpi(le) = rpi(le)**rk                                                     
 30   continue                                                                  
c                                                                               
      do 6 le=1,kdimpi                                                          
      print 100, le, ci(le), si(le)                                             
100   format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
    6 continue                                                                  
c                                                                               
      print 200                                                                 
  200 format (1h0)                                                              
      do 7 le=1,kdimi                                                           
      print 300, le, cl(le), sl(le), del(le)                                    
  300 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,          
     1 'del=', f6.3)                                                            
    7 continue                                                                  
c                                                                               
      print 400, (rpi(le), le=1,kdimmi)                                         
  400 format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
c                                                                               
      return                                                                    
      end                                                                       
