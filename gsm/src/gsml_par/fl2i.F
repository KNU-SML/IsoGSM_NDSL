#define JCAP1 jcap1
#define JCAP jcap
#define JCAPP jcapp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT2 lnt2
#define LNT lnt
      subroutine fl2i(fp,fm,fln,qln,llstr,llens,lwvdef,lota)
#include <paramodel.h>
c 
c  yifeng cui, off rsm mp   02-09-05
c
#ifdef RSM
#undef MP
#endif
#ifdef MP
#define JCAPS jcapp_
#define LLN2S lln2p_
#define LLN22S lln22p_
#else
#define JCAPS jcap_
#define LLN2S lnt2_
#define LLN22S lnt22_
#endif
      parameter (len0p=JCAPS)                                                   
      parameter (len0m=JCAPS)                                                   
      parameter (lnt=lnt_)                                                      
      parameter (jcap=jcap_)                                                    
      dimension fp(2,0:len0p,lota), fm(2,0:len0m,lota),                         
     *          qln(LLN2S), fln(LLN22S,lota)
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
      offset(n,l) = (jcap+1)*(jcap+2) - (jcap+1-l)*(jcap+2-l) + 2*(n-l)         
c                                                                               
c     ----------------------------------------------------------------          
c     compute the coefficients of the expansion in spherical harmonics          
c     of the field at each level                                                
c     ----------------------------------------------------------------          
c                                                                               
      lntx=0
      do 1000 lx=1,llens
        l=lwvdef(llstr+lx)
        lnt0=offset(l,l)-lntx
        lntx=lntx+offset(l+1,l+1)-offset(l,l)
        ll=lx-1
c
#define DEFAULT
#ifdef SGERX1
#undef DEFAULT
         ls=l*((2*jcap+3)-l)                                                    
         lls=ls-lnt0
#endif
c                                                                               
c        compute the even (n-l) expansion coefficients for each level           
c        ------------------------------------------------------------           
c                                                                               
c        real part                                                              
c                                                                               
#ifdef SGERX1
         call sgerx1((jcap+2-l)/2,lota,1.,qln(lls+1),4,                          
     1               fp(1,ll,1),(len0p+1)*2,fln(lls+1,1),4,LLN22S)
#endif
#ifdef DEFAULT
         do n = l, jcap, 2                                                      
            do k = 1, lota                                                      
               nl=offset(n,l)-lnt0
               fln(nl+1,k) = fln(nl+1,k)
     1                      + fp(1,ll,k)*qln(nl+1)
            end do                                                              
         end do                                                                 
#endif
c                                                                               
c        imaginary part                                                         
c                                                                               
#ifdef SGERX1
         call sgerx1((jcap+2-l)/2,lota,1.,qln(lls+2),4,                          
     1               fp(2,ll,1),(len0p+1)*2,fln(lls+2,1),4,LLN22S) 
#endif
#ifdef DEFAULT
         do n = l, jcap, 2                                                      
            do k = 1, lota                                                      
               nl=offset(n,l)-lnt0
               fln(nl+2,k) = fln(nl+2,k)
     1                      + fp(2,ll,k)*qln(nl+2)
            end do                                                              
         end do                                                                 
#endif
c                                                                               
c        compute the odd (n-l) expansion coefficients for each level            
c        -----------------------------------------------------------            
#ifdef SGERX1
         if(l.lt.jcap) then                                                     
#endif
c                                                                               
c        real part                                                              
c                                                                               
#ifdef SGERX1
           call sgerx1((jcap+1-l)/2,lota,1.,qln(lls+3),4,                        
     1                 fm(1,ll,1),(len0m+1)*2,fln(lls+3,1),4,LLN22S)               
#endif
#ifdef DEFAULT
           do n = l+1, jcap, 2                                                  
              do k = 1, lota                                                    
                nl=offset(n,l)-lnt0
                fln(nl+1,k) = fln(nl+1,k)
     1                      + fm(1,ll,k)*qln(nl+1)          
              end do                                                            
           end do                                                               
#endif
c                                                                               
c        imaginary part                                                         
c                                                                               
#ifdef SGERX1
           call sgerx1((jcap+1-l)/2,lota,1.,qln(lls+4),4,                        
     1                 fm(2,ll,1),(len0m+1)*2,fln(lls+4,1),4,LLN22S)               
#endif
#ifdef DEFAULT
           do n = l+1, jcap, 2                                                  
              do k = 1, lota                                                    
                nl=offset(n,l)-lnt0
                fln(nl+2,k) = fln(nl+2,k)
     1                       + fm(2,ll,k)*qln(nl+2)
              end do                                                            
           end do                                                               
#endif
c                                                                               
#ifdef SGERX1
         endif                                                                  
#endif
c
 1000 continue
c                                                                               
      return                                                                    
      end                                                                       
