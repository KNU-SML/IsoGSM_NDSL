#define JCAP1 jcap1
#define JCAP jcap
#define JCAPP jcapp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT2 lnt2
#define LNT lnt
      subroutine fl2ip(fp,fm,fln,qln,llstr,llens,lwvdef,lota)
c                                                                               
c parallelize at k level.                                                       
c                                                                               
#include <paramodel.h>
#ifdef MP
#define JCAPS jcapp_
#define LLN2S lln2p_
#define LLN22S lln22p_
#else
#define JCAPS jcap_
#define LLN2S lnt2_
#define LLN22S lnt22_
#endif
      parameter (len0p=JCAPS )
      parameter (len0m=JCAPS )
      parameter (lnt=lnt_ )                                                     
      parameter (jcap=jcap_ )                                                   
      dimension fp(2,0:len0p,lota), fm(2,0:len0m,lota),                         
     .          qln(LLN2S), fln(LLN22S,lota)
      dimension lwvdef(jcap1_)
c                                                                               
c     local scalars                                                             
c     -------------                                                             
c                                                                               
      integer n, l, k                                                           
c                                                                               
c     statement function                                                        
c     ------------------                                                        
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
c               |   x x x x                                                     
c               v   x x x x x                                                   
c                   x x x x x x                                                 
c                                                                               
c          order of the matrix elements in memory:                              
c                                                                               
c          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
c                                                                               
      integer offset                                                            
      offset(n,l)=(jcap+1)*(jcap+2)-(jcap+1-l)*(jcap+2-l)+2*(n-l)
c                                                                               
c     ----------------------------------------------------------------          
c     compute the coefficients of the expansion in spherical harmonics          
c     of the field at each level                                                
c     ----------------------------------------------------------------          
c                                                                               
c
#ifdef CRAY_THREAD
ccfpp$ cncall                                                                   
#endif
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(fp,fm,fln,qln,llstr,llens,lwvdef,jcap,lota)
c$& local(k,l,n,ll,nl,lntx,lnt0,lx)                                                                
#endif
#ifdef OPENMP
!$omp parallel do private(k,l,n,ll,nl,lntx,lnt0,lx)
#endif
#ifdef CRAY_THREAD
cmic$ doall shared(fp,fm,fln,qln,llstr,llens,lwvdef,jcap,lota)
cmic$1      private(k,l,n,ll,nl,lntx,lnt0,lx)
#endif
c                                                                               
      do k=1,lota  ! mj 5/8/1998                                                
c
        lntx=0
        do lx=1,llens
          l=lwvdef(llstr+lx)
          lnt0=offset(l,l)-lntx
          lntx=lntx+offset(l+1,l+1)-offset(l,l)
          ll=lx-1
c                                                                               
c        compute the even (n-l) expansion coefficients for each level           
c        ------------------------------------------------------------           
c                                                                               
c        real part                                                              
c                                                                               
          do n = l, jcap, 2                                                      
            nl=offset(n,l)-lnt0
            fln(nl+1,k) = fln(nl+1,k)
     1                   + fp(1,ll,k)*qln(nl+1)
          end do                                                                 
c                                                                               
c        imaginary part                                                         
c                                                                               
          do n = l, jcap, 2                                                      
            nl=offset(n,l)-lnt0
            fln(nl+2,k) = fln(nl+2,k)
     1                   + fp(2,ll,k)*qln(nl+2)
          end do                                                                 
c                                                                               
c        compute the odd (n-l) expansion coefficients for each level            
c        -----------------------------------------------------------            
c                                                                               
c        real part                                                              
c                                                                               
          do n = l+1, jcap, 2                                                  
            nl=offset(n,l)-lnt0
            fln(nl+1,k) = fln(nl+1,k)
     1                   + fm(1,ll,k)*qln(nl+1)
          end do                                                               
c                                                                               
c        imaginary part                                                         
c                                                                               
          do n = l+1, jcap, 2                                                  
            nl=offset(n,l)-lnt0
            fln(nl+2,k) = fln(nl+2,k)
     1                   + fm(2,ll,k)*qln(nl+2)
          end do                                                               
c                                                                               
        enddo                                                                   
      enddo     !mj 5/8/1998                                                    
c                                                                               
      return                                                                    
      end                                                                       
