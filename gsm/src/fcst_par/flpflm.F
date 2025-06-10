#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
#define LCAP22P lcap22p
#define LCAPP lcapp
#define LONF22 lonf22
#define LONF lonf
      subroutine flpflm(flp,flm,anl,lwvlen,lota)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c.................................................................              
c  for loopa                                                                    
c                                                                               
#include <paramodel.h>
#ifdef MP
#define JCAP1S jcap1p_
#define LCAPS  lcapp_
#define LCAP22S lcap22p_
#else
#define JCAP1S jcap1_
#define LCAPS lonf_
#define LCAP22S lonf22_
#endif
       dimension                                                                
     x flp(2,JCAP1S,lota),flm(2,JCAP1S,lota),                             
     x anl(LCAP22S,lota)                                                     
c                                                                               
c.................................................................              
c                                                                               
      do k=1,lota                                                            
      do  ll=1,lwvlen
c do n.hemi                                                                     
       flp(1,ll,k)=anl(2*(ll-1)+1,k)+                                           
     1             anl(2*(ll-1)+1+LCAPS,k)
       flp(2,ll,k)=anl(2*(ll-1)+2,k)+                                           
     1             anl(2*(ll-1)+2+LCAPS,k)
c                                                                               
c do s.hemi                                                                     
c                                                                              
       flm(1,ll,k)=anl(2*(ll-1)+1,k)-                                           
     1             anl(2*(ll-1)+1+LCAPS,k)                                      
       flm(2,ll,k)=anl(2*(ll-1)+2,k)-                                           
     1             anl(2*(ll-1)+2+LCAPS,k)
      enddo                                                                     
      enddo                                                                     
c.................................................................              
      return                                                                    
      end                                                                       
