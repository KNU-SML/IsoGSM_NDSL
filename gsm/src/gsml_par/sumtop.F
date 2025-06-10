#define JCAP1 jcap1
#define JCAP jcap
#define JCAPP jcapp
#define LCAPP lcapp
#define LEVS levs
#define LEVSP levsp
#define LONF lonf
      subroutine sumtop(ap,top,qvv,llstr,llens,lwvdef,kall) 
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
cc                                                                              
#include <paramodel.h>
#ifdef RSM
#undef MP
#endif
#ifdef MP
#define JCAPS jcapp_
#define LEVSS levsp_
#define LCAPS lcapp_
#else
#define JCAPS jcap_
#define LEVSS levs_
#define LCAPS lonf_
#endif
      parameter (len0=LCAPS)
      parameter (lenh=len0/2)
      parameter (lev2p1=2*LEVSS+1)                                          
cc                                                                              
      dimension     ap(2,0:len0,kall)                                           
      dimension    top(2,0:JCAPS,kall)                                          
      dimension    qvv(2,0:JCAPS)                                                  
      dimension    lwvdef(jcap1_)                                                  
cc local array
      dimension   qtop(2,0:JCAPS)
c
cc                                                                              
      do 1000 lx=1,llens
        l=lwvdef(llstr+lx)
        ll=lx-1
c
        do 500 k=1,kall                                                           
c
          qtop(1,ll) = top(1,ll,k) * qvv(1,ll)
          qtop(2,ll) = top(2,ll,k) * qvv(2,ll)
c
          ap(1,ll,k) = ap(1,ll,k) + qtop(1,ll)
          ap(2,ll,k) = ap(2,ll,k) + qtop(2,ll)
cc                                                                              
          if( mod(l,2).eq.0 ) then
c
c odd wavenumber  for l=0,2,4,.....
c
            ap(1,ll+lenh,k) = ap(1,ll+lenh,k) - qtop(1,ll)
            ap(2,ll+lenh,k) = ap(2,ll+lenh,k) - qtop(2,ll)
cc  
          else
cc                                                                              
c even wavenumber for l=1,3,5,......
c
            ap(1,ll+lenh,k) = ap(1,ll+lenh,k) + qtop(1,ll)
            ap(2,ll+lenh,k) = ap(2,ll+lenh,k) + qtop(2,ll)
c
          endif
cc                                                                              
  500   continue                                                                  
c
 1000 continue
cc                                                                              
      return                                                                    
      end                                                                       
