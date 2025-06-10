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
      subroutine newsig(ci, si, del, sl, cl, rpi)                               
#include <paramodel.h>
      save                                                                      
c ------------------********************-----------------------                 
c this routine sets coordinates for levels and use phillips method              
c to get layers. ( can be stretched levels or not )                             
c  --------------------------------------------------------------               
c                                                                               
#include <paramter.h>
c                                                                               
      parameter(kdimp=kdim+1,kdimm=kdim-1)                                      
c                                                                               
      dimension ci(kdimp),si(kdimp),del(kdim),sl(kdim),cl(kdim),                
     &          rpi(kdimm)                                                      
c                                                                               
      real rk,rk1,rkinv                                                         
c                                                                               
#include <modlsigs.h>
c                                                                               
#include <parmchgr.h>
      parameter(kdimpi=kdimi+1,kdimmi=kdimi-1)                                  
      common/inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),                     
     1               ciin(kdimpi),clin(kdimi),rpiin(kdimmi)                     
     2               ,ak(kdimpi),bk(kdimpi)                                     ! 20070501(12)

c                                                                               
      data ifp/0/                                                               
c***********************************************************************        
c                                                                               
      ci(1) = 0.                                                                
      do 2 k=1,kdim                                                             
      del(k)=delmdl(k)                                                          
      ci(k+1)=ci(k)+delmdl(k)                                                   
    2 continue                                                                  
      ci(kdimp)=1.                                                              
c                                                                               
      rk  = rd_/cp_                                                             
      rk1 = rk + 1.                                                             
      rkinv=1./rk                                                               
c                                                                               
      levs=kdim                                                                 
c                                                                               
      do li=1,kdimp                                                             
      si(li) = 1. - ci(li)                                                      
      enddo                                                                     
c                                                                               
      do 4 le=1,kdim                                                            
      dif = si(le)**rk1 - si(le+1)**rk1                                         
      dif = dif / (rk1*(si(le)-si(le+1)))                                       
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
