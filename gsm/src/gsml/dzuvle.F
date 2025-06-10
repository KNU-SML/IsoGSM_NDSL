#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
#define LEVS levs
#define LEVSP levsp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT2 lnt2
#define RERTH rerth
      subroutine dzuvle(d,z,u,v,utop,vtop,llstr,llens,lwvdef)
c
c   yifeng cui, off rsm mp,   02-09-05
c      
#include <paramodel.h>
      parameter (jcap=jcap_)                                                    
      parameter (lnep=(jcap+2)*(jcap+3)/2)                                      
c
#ifdef RSM
#undef MP
#endif
#ifdef MP
#define JCAP1S jcap1p_
#define LEVSS levsp_
#define LLN2S lln2p_
#define LLN22S lln22p_
#else
#define JCAP1S jcap1_
#define LEVSS levs_
#define LLN2S lnt2_
#define LLN22S lnt22_
#endif
      dimension  d(LLN22S,LEVSS),z(LLN22S,LEVSS)
      dimension  u(LLN22S,LEVSS),v(LLN22S,LEVSS)
      dimension  utop(2,JCAP1S,LEVSS)
      dimension  vtop(2,JCAP1S,LEVSS)
      dimension  lwvdef(jcap1_)                                                        
      dimension  e(lnep)                                                        
c                                                                               
      save e, ifirst                                                            
c                                                                               
c   array e =eps/n                                                              
c   array e =eps/n                                                              
c            eps/n=0. for n=l                                                   
c   array e =eps/n                                                              
c   array e =eps/n                                                              
c                                                                               
      je(n,l) =((jcap+2)*(jcap+3)-(jcap+2-l)*(jcap+3-l))/2+n-l                  
c                                                                               
      jc(n,l) = (jcap+1)*(jcap+2)-(jcap+1-l)*(jcap+2-l)+2*(n-l)                 
c                                                                               
      data ifirst/1/                                                            
      if(ifirst.ne.1)go to 280                                                  
c ------------------ initial setting --------------------                       
      do 220 l=0,jcap                                                           
             n=l                                                                
      e(je(n,l)+1)=0.                                                           
220   continue                                                                  
      do 260 l=  0,jcap                                                         
      do 240 n=l+1,jcap+1                                                       
      rn=n                                                                      
      rl=l                                                                      
      a=(rn*rn-rl*rl)/(4.*rn*rn-1.)                                             
      e(je(n,l)+1)=sqrt(a) / rn                                                 
240   continue                                                                  
260   continue                                                                  
      ifirst=0                                                                  
280   continue                                                                  
c -------------------- end of initial setting ---------------
      lntx=0
      do 2000 ll=1,llens
        l=lwvdef(llstr+ll)
        lnt0=jc(l,l)-lntx
        lntx=lntx+jc(l+1,l+1)-jc(l,l)
c
#ifdef ORIGIN_THREAD
c$doacross share(d,z,u,v,utop,vtop,e,lnt0,l,ll),
c$&        local(k,n,rl,rn,j,jc0,je0)
#endif
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 shared(d,z,u,v,utop,vtop,e,lnt0,l,ll)
cmic$1 private(k,n,rl,rn,j,jc0,je0)
#endif
#ifdef OPENMP
!$omp parallel do private(k,n,rl,rn,j,jc0,je0)
#endif
      do 1000 k=1,LEVSS                                                        
c
      if(l.eq.0) then
        do 320 n=0,jcap
          jc0=jc(n,l)
          u(jc0+1,k)=0.0                                                        
          u(jc0+2,k)=0.0                                                        
c                                                                               
          v(jc0+1,k)=0.0 
          v(jc0+2,k)=0.0
320     continue  
      endif
c                                                                               
c----- l=1,jcap
c 
      if( l.ge.1 ) then
          rl=l                                                                  
        do 420 n=l,jcap                                                           
          rn=n                                                                  
          jc0=jc(n,l)-lnt0
c                                                                               
          u(jc0+2,k)=-rl*d(jc0+1,k)/(rn*(rn+1.))                            
          u(jc0+1,k)= rl*d(jc0+2,k)/(rn*(rn+1.))                            
c                                                                               
          v(jc0+2,k)=-rl*z(jc0+1,k)/(rn*(rn+1.))                            
          v(jc0+1,k)= rl*z(jc0+2,k)/(rn*(rn+1.))                            
420     continue                                                                  
      endif
c                                                                               
c----- l=  0,jcap-1
c 
      if( l.le.jcap-1 ) then
        do 520 n=l+1,jcap                                                         
          jc0=jc(n,l)-lnt0
          je0=je(n,l)
c
          u(jc0+1,k)=u(jc0+1,k)-e(je0+1)*z(jc0-1,k)                 
          u(jc0+2,k)=u(jc0+2,k)-e(je0+1)*z(jc0  ,k)                 
c                                                                               
          v(jc0+1,k)=v(jc0+1,k)+e(je0+1)*d(jc0-1,k)                 
          v(jc0+2,k)=v(jc0+2,k)+e(je0+1)*d(jc0  ,k)                 
520     continue                                                                  
        do 620 n=l,jcap-1                                                         
          jc0=jc(n,l)-lnt0
          je0=je(n+1,l)
c
          u(jc0+1,k)=u(jc0+1,k)+e(je0+1)*z(jc0+3,k)               
          u(jc0+2,k)=u(jc0+2,k)+e(je0+1)*z(jc0+4,k)               
c                                                                               
          v(jc0+1,k)=v(jc0+1,k)-e(je0+1)*d(jc0+3,k)               
          v(jc0+2,k)=v(jc0+2,k)-e(je0+1)*d(jc0+4,k)               
620     continue                                                                  
      endif
c                                                                               
      n=jcap+1                                                                  
c
c----- l=0,jcap
c     
      jc0=jc(n,l)-lnt0
      je0=je(n,l)
c
      utop(1,ll,k)=-e(je0+1)*z(jc0-1,k)                                
      utop(2,ll,k)=-e(je0+1)*z(jc0  ,k)                                
c                                                                               
      vtop(1,ll,k)= e(je0+1)*d(jc0-1,k)                                
      vtop(2,ll,k)= e(je0+1)*d(jc0  ,k)                                
c                                                                               
c----- l=0,jcap_
c  
      utop(1,ll,k)=utop(1,ll,k)*rerth_                                            
      utop(2,ll,k)=utop(2,ll,k)*rerth_                                            
      vtop(1,ll,k)=vtop(1,ll,k)*rerth_                                            
      vtop(2,ll,k)=vtop(2,ll,k)*rerth_                                            
c
1000  continue                                                                  
c
2000  continue
c
      do k=1,LEVSS
        do j=1,LLN2S
          u(j,k)=u(j,k)*rerth_ 
          v(j,k)=v(j,k)*rerth_
        enddo
      enddo
c                                                                               
      return                                                                    
      end                                                                       
