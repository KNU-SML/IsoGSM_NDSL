#define JCAP1 jcap1
#define JCAP jcap
#define JCAPP jcapp
#define LCAPP lcapp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT2 lnt2
#define LNT lnt
#define LONF lonf
      subroutine sums2i(fln,ap,qln,llstr,llens,lwvdef,levs)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
c 
c if rsm, run se
c
#ifdef RSM
#undef MP
#endif
#ifdef MP
#define LLN2S lln2p_
#define LLN22S lln22p_
#define LCAPS lcapp_
#else
#define LLN2S lnt2_
#define LLN22S lnt22_
#define LCAPS lonf_
#endif
      parameter (len0=LCAPS)
      parameter (lenh=len0/2)                                                  
      parameter (lnt=lnt_)                                                      
      parameter (jcap=jcap_)                                                    
      dimension ap(2,0:len0,levs), qln(LLN2S), fln(LLN22S,levs)
      dimension lwvdef(jcap1_)
c                                                                               
c     local scalars                                                             
c     -------------                                                             
c                                                                               
      integer i, n, l, k                                                        
      real    evenr, eveni                                                      
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
c               |   x x x x                                                     
c               v   x x x x x                                                   
c                   x x x x x x                                                 
c                                                                               
c          order of the matrix elements in memory:                              
c                                                                               
c          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
c                                                                               
      integer offset                                                            
      offset(n,l)=(jcap+1)*(jcap+2)-(jcap-l+1)*(jcap-l+2)+2*(n-l)
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
      term(i,n,l,k) = qln(offset(n,l)-lnt0+i)*fln(offset(n,l)-lnt0+i,k)
c                                                                               
c     zero the accumulators                                                     
c     ---------------------                                                     
c                                                                               
c     lens=l1e-l1s+l2e-l2s+1
      do k = 1, levs                                                            
         do l = 0, len0
            ap(1,l,k) = 0.                                                      
            ap(2,l,k) = 0.                                                      
         end do                                                                 
      end do                                                                    
c                                                                               
c     compute the even and odd (n-l) components                                 
c     of the fourier coefficients                                               
c     ---------------------------------------------------------                 
c                                                                               
      lntx=0
      do 1000 lx=1,llens
        l=lwvdef(llstr+lx)
        lnt0=offset(l,l)-lntx
        lntx=lntx+offset(l+1,l+1)-offset(l,l)
        ll=lx-1
c
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
         ls=l*((2*jcap+3)-l)                                                    
         lls=ls-lnt0
#endif
c                                                                               
c        compute the sum of the even (n-l) terms for each level                 
c        ------------------------------------------------------                 
c                                                                               
c        real part                                                              
c                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
         call sgemvx1(levs,(jcap+2-l)/2,1.,fln(lls+1,1),LLN22S,4,
     1                qln(lls+1),4,1.,ap(1,ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
         do n = l, jcap, 2                                                      
            do k = 1, levs
c               print *, 'n=',n,'l=',l,'k=',k,'term(1,n,l,k)=',term(1,n,l,k),'ap(1,ll,k)=',ap(1,ll,k)                                                      
               ap(1,ll,k) = ap(1,ll,k) + term(1,n,l,k)
            end do                                                              
         end do                                                                 
#endif
c                                                                               
c        imaginary part                                                         
c                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
         call sgemvx1(levs,(jcap+2-l)/2,1.,fln(lls+2,1),LLN22S,4,
     1                qln(lls+2),4,1.,ap(2,ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
         do n = l, jcap, 2                                                      
            do k = 1, levs                                                      
               ap(2,ll,k) = ap(2,ll,k) + term(2,n,l,k)
            end do                                                              
         end do                                                                 
#endif
c                                                                               
c        compute the sum of the odd (n-l) terms for each level                  
c        -----------------------------------------------------                  
c                                                                               
c          real part                                                            
c                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
         if(l.lt.jcap) then                                                     
           call sgemvx1(levs,(jcap+1-l)/2,1.,fln(lls+3,1),LLN22S,4,
     1                qln(lls+3),4,1.,ap(1,lenh+ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
           do n = l+1, jcap, 2                                                  
              do k = 1, levs                                                    
                 ap(1,lenh+ll,k) = ap(1,lenh+ll,k) + term(1,n,l,k)
              end do                                                            
           end do                                                               
#endif
c                                                                               
c          imaginary part                                                       
c                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
           call sgemvx1(levs,(jcap+1-l)/2,1.,fln(lls+4,1),LLN22S,4,
     1                qln(lls+4),4,1.,ap(2,lenh+ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
           do n = l+1, jcap, 2                                                  
              do k = 1, levs                                                    
                 ap(2,lenh+ll,k) = ap(2,lenh+ll,k) + term(2,n,l,k)
              end do                                                            
           end do                                                               
#endif
c                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
         endif                                                                  
#endif
c
c     compute the fourier coefficients for each level                           
c     -----------------------------------------------                           
c                                                                               
        do k = 1, levs                                                            
          evenr = ap(1,ll,k)                                                   
          eveni = ap(2,ll,k)                                                   
          ap(1,ll,k) = ap(1,ll,k) + ap(1,lenh+ll,k)                              
          ap(2,ll,k) = ap(2,ll,k) + ap(2,lenh+ll,k)                              
          ap(1,lenh+ll,k) = evenr - ap(1,lenh+ll,k)                             
          ap(2,lenh+ll,k) = eveni - ap(2,lenh+ll,k)                             
        end do                                                                    
c
 1000 continue
c                                                                               
      return                                                                    
      end                                                                       
