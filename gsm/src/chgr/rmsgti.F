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
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine rmsgti(gz,q,x,y,w,del,r)                                          
#include <paramodel.h>
      save                                                                      
c                                                                               
#include <parmchgr.h>
#include <paramter.h>
c                                                                               
      parameter(mdim=(mwave+1)*(mwave+2))                                       
c                                                                               
      dimension gz(mdim),q(mdim),x(mdim,kdimi),y(mdim,kdimi)                             
      dimension w(mdim,kdimi),r(mdim,kdimi)                                     
      dimension del(kdimi)                                                      
c                                                                               
      dimension rx(kdimi),ry(kdimi),rw(kdimi),rr(kdimi)                         
c                                                                               
      vr=0.                                                                     
      vx=0.                                                                     
      vy=0.                                                                     
      vw=0.                                                                     
      do 10 k=1,kdimi                                                           
      call hbartr(x(1,k),x(1,k),rx(k))                                          
      call hbartr(y(1,k),y(1,k),ry(k))                                          
      call hbartr(w(1,k),w(1,k),rw(k))                                          
      call hbartr(r(1,k),r(1,k),rr(k))                                          
      vx=vx+rx(k)*del(k)                                                        
      vy=vy+ry(k)*del(k)                                                        
      vw=vw+rw(k)*del(k)                                                        
      vr=vr+rr(k)*del(k)                                                        
10    continue                                                                  
c                                                                               
      call hbartr(q,q,rq)                                                       
      call hbartr(gz,gz,rgz)                                                       
      print 100,vx,vw,vy,vr,rq,rgz                                                  
100   format('rmsgti div vort temp mixratio ln(ps) orog',6(e8.3,1x))                 
200   format(1h ,4(2x,e8.3))                                                    
      do 40 k=1,kdimi                                                           
      print 200,rx(k),rw(k),ry(k),rr(k)                                         
40    continue                                                                  
      return                                                                    
      end                                                                       
