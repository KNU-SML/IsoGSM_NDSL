#define JCAP1 jcap1
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
      subroutine uvtodz(uln,vln,di,ze,topuln,topvln,llstr,llens,lwvdef)
#include <paramodel.h>
#ifdef MP
#define JCAPS jcapp_
#define LLN2S lln2p_
#define LLN22S lln22p_
#define LEVSS levsp_
#else
#define JCAPS jcap_
#define LLN2S lnt2_
#define LLN22S lnt22_
#define LEVSS levs_
#endif
      parameter (jcap=jcap_)                                                    
      dimension  di(LLN22S,LEVSS)
      dimension  ze(LLN22S,LEVSS)
      dimension uln(LLN22S,LEVSS)
      dimension vln(LLN22S,LEVSS)
      dimension eps(lnt22_)      
      dimension   topuln(2,0:JCAPS,LEVSS)
      dimension   topvln(2,0:JCAPS,LEVSS)
      dimension   topeps(0:jcap_)        
      dimension   lwvdef(jcap1_)        
c                                                                               
      save ifirst,eps,topeps                                                    
c                                                                               
c     local scalars                                                             
c     -------------                                                             
c                                                                               
      integer i, n, l, k                                                        
c                                                                               
c     statement functions                                                       
c     -------------------                                                       
c                                                                               
c     offset(n,l) is the offset in words                                        
c     to the (n,l)-element of a lower                                           
c     triangular matrix of complex numbers                                      
c     in an array containing the matrix                                         
c     packed in column-major order,                                             
c     where l and n range from 0 to jcap,                                       
c     inclusive                                                                 
c                                                                               
c          lower triangular matrix of complex numbers:                          
c                                                                               
c                     l -->                                                     
c                                                                               
c                   x                                                           
c               n   x x                                                         
c                   x x x                                                       
c               ?   x x x x                                                     
c               v   x x x x x                                                   
c                   x x x x x x                                                 
c                                                                               
c          order of the matrix elements in memory:                              
c                                                                               
c          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
c                                                                               
      integer offset                                                            
      offset(n,l) = (jcap+1)*(jcap+2) - (jcap-l+1)*(jcap-l+2) + 2*(n-l)         
c                                                                               
c     ---                                                                       
c                                                                               
c     term(1,n,l,k) and term(2,n,l,k) are                                       
c     the real and imaginary part, resp.,                                       
c     of exp((0,1)*l*phi) times the (n,l) term                                  
c     in the expansion in spherical                                             
c     harmonics of the field at level k,                                        
c     where phi is the azimuthal angle                                          
c                                                                               
c     term(i,n,l,k) = di(offset(n,l)+i,k)                                       
c                                                                               
c                                                                               
      data ifirst/0/                                                            
      if(ifirst.ne.0) go to 999                                                 
      do l = 0, jcap                                                            
         do n = l, jcap                                                         
            temp=((n*n-l*l)/(4.*n*n-1.))                                        
            if(n.eq.0) temp=0.                                                  
            temp=sqrt(temp)                                                     
            eps(offset(n,l)+1)=temp                                             
            eps(offset(n,l)+2)=temp                                             
         end do                                                                 
            n=jcap+1                                                            
            temp=((n*n-l*l)/(4.*n*n-1.))                                        
            topeps(l)=sqrt(temp)                                                
      end do                                                                    
c                                                                               
      ifirst=1                                                                  
999   continue                                                                  
c --------------------------------------------------------------
c
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$1 autoscope                                                                
#endif
      do 1000 k = 1, LEVSS                                        
c
        lntx=0
        do lx=1,llens
          l=lwvdef(llstr+lx)
          lnt0=offset(l,l)-lntx
          lntx=lntx+offset(l+1,l+1)-offset(l,l)
          ll=lx-1
c
c        the case n=l                                                           
c-----do l = 1, jcap-1                                                          
c 
          if( l.ge.1 .and. l.le.jcap-1 ) then
            rl= l
            n = l                                                               
            rn= n                                                                  
            nl=offset(n,l)-lnt0
            npl=offset(n+1,l)-lnt0
c
            ze(nl+1,k)=-rl*vln(nl+2,k)                       
     1      -rn*eps(offset(n+1,l)+1)*uln(npl+1,k)                    
c                                                                               
            ze(nl+2,k)= rl*vln(nl+1,k)                       
     1      -rn*eps(offset(n+1,l)+2)*uln(npl+2,k)                    
c                                                                               
            di(nl+1,k)=-rl*uln(nl+2,k)                       
     1      +rn*eps(offset(n+1,l)+1)*vln(npl+1,k)                    
c                                                                               
            di(nl+2,k)= rl*uln(nl+1,k)                       
     1      +rn*eps(offset(n+1,l)+2)*vln(npl+2,k)                    
c                                                                               
          endif
c 
c---  do l=0
c 
          if( l.eq.0 ) then
            ze(1,k)=0.                                                               
            ze(2,k)=0.                                                               
            di(1,k)=0.                                                               
            di(2,k)=0.                                                               
          endif
c
c---  do l = 0, jcap                                                            
c  
          rl=l
          do n = l+1, jcap-1                                                     
            rn=n                                                                  
            nl=offset(n,l)-lnt0
            npl=offset(n+1,l)-lnt0
            nml=offset(n-1,l)-lnt0
c        
            ze(nl+1,k)=-rl*vln(nl+2,k)                       
     1            -rn*eps(offset(n+1,l)+1)*uln(npl+1,k)                    
     2       +(rn+1.)*eps(offset(n  ,l)+1)*uln(nml+1,k)                 
c                                                                               
            ze(nl+2,k)= rl*vln(nl+1,k)                       
     1            -rn*eps(offset(n+1,l)+2)*uln(npl+2,k)                    
     2       +(rn+1.)*eps(offset(n  ,l)+2)*uln(nml+2,k)                 
c                                                                               
            di(nl+1,k)=-rl*uln(nl+2,k)                       
     1            +rn*eps(offset(n+1,l)+1)*vln(npl+1,k)                    
     2       -(rn+1.)*eps(offset(n  ,l)+1)*vln(nml+1,k)                 
c                                                                               
            di(nl+2,k)= rl*uln(nl+1,k)                       
     1            +rn*eps(offset(n+1,l)+2)*vln(npl+2,k)                    
     2       -(rn+1.)*eps(offset(n  ,l)+2)*vln(nml+2,k)                 
c                                                                               
          end do                                                                 
c   
c do top row involving u,v at n=jcap+1                                          
ccccccc                                                                         
          n =  jcap                                                             
          rn=n                                                                  
c
c-----do l = 0, jcap                                                            
c
          rl=l                                                                     
          nl=offset(n,l)-lnt0
          nml=offset(n-1,l)-lnt0
c     
          ze(nl+1,k)=-rl*vln(nl+2,k)                       
     2    +(rn+1.)*eps(offset(n,l)+1)*uln(nml+1,k)                 
c                                                                               
          ze(nl+2,k)= rl*vln(nl+1,k)                       
     2    +(rn+1.)*eps(offset(n,l)+2)*uln(nml+2,k)                 
c                                                                               
          di(nl+1,k)=-rl*uln(nl+2,k)                       
     2    -(rn+1.)*eps(offset(n,l)+1)*vln(nml+1,k)                 
c                                                                               
          di(nl+2,k)= rl*uln(nl+1,k)                       
     2    -(rn+1.)*eps(offset(n,l)+2)*vln(nml+2,k)                 
c                                                                               
ccccccc                                                                         
c-----do l = 0, jcap                                                            
c  
c         nl=offset(n,l)-lnt0
c
          ze(nl+1,k)=ze(nl+1,k)-rn*topeps(l)*topuln(1,ll,k)
c                                                                               
          ze(nl+2,k)=ze(nl+2,k)-rn*topeps(l)*topuln(2,ll,k)
c                                                                               
          di(nl+1,k)=di(nl+1,k)+rn*topeps(l)*topvln(1,ll,k)
c                                                                               
          di(nl+2,k)=di(nl+2,k)+rn*topeps(l)*topvln(2,ll,k)
c                                                                               
        enddo                                                                    
ccccccc                                                                         
1000  continue                                                                  
c
      do k=1,LEVSS
        do j=1,LLN2S
          di(j,k)=di(j,k)/rerth_                                            
          ze(j,k)=ze(j,k)/rerth_                                            
        end do                                                             
      enddo
c 
      return                                                                    
      end                                                                       
