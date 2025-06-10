#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NTRAC ntrac

      subroutine lfmfilt(ifstep,thour)                                          
c                                                                               
#include <paramodel.h>
#include <comfspec.vr.h>
#include <comfphys.h>
#include <comlfm.h>
c                                                                               
#ifndef NOPRINT
      write(6,*) 'doing lfmfilt for ifstep=',ifstep                             
#endif
      do l=1,lnt2_                                                              
        fq(l)=fq(l)+q(l)*weight(ifstep)                                         
      enddo                                                                     
      do k=1,levs_                                                              
        do l=1,lnt2_                                                            
          fte(l,k)=fte(l,k)+te(l,k)*weight(ifstep)                              
          fdi(l,k)=fdi(l,k)+di(l,k)*weight(ifstep)                              
          fze(l,k)=fze(l,k)+ze(l,k)*weight(ifstep)                              
        enddo                                                                   
      enddo                                                                     
      do k=1,levh_                                                              
        do l=1,lnt2_                                                            
          frq(l,k)=frq(l,k)+rq(l,k)*weight(ifstep)                              
        enddo                                                                   
      enddo                                                                     
c                                                                               
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          ftsea(i,j)=ftsea(i,j)+tsea(i,j)*weight(ifstep)                        
          fsheleg(i,j)=fsheleg(i,j)+sheleg(i,j)*weight(ifstep)                  
          ftg3(i,j)=ftg3(i,j)+tg3(i,j)*weight(ifstep)                           
          fzorl(i,j)=fzorl(i,j)+zorl(i,j)*weight(ifstep)                        
          fplantr(i,j)=fplantr(i,j)+plantr(i,j)*weight(ifstep)                  
          fcv(i,j)=fcv(i,j)+cv(i,j)*weight(ifstep)                              
          do il = 1, 4
          falbedo(i,j,il)=falbedo(i,j,il)+albedo(i,j,il)*weight(ifstep)                  
          enddo
          ff10m(i,j)=ff10m(i,j)+f10m(i,j)*weight(ifstep)                        
          fcanopy(i,j)=fcanopy(i,j)+canopy(i,j)*weight(ifstep)                  
          isl=nint(slmsk(i,j))+1                                                
          islmsk(i,j,isl)=islmsk(i,j,isl)+1                                     
          if(cvb(i,j).ne.cvb0) then                                             
            fcvb(i,j)=fcvb(i,j)+cvb(i,j)*weight(ifstep)                         
            wcvb(i,j)=wcvb(i,j)+weight(ifstep)                                  
          endif                                                                 
          if(cvt(i,j).ne.cvt0) then                                             
            fcvt(i,j)=fcvt(i,j)+cvt(i,j)*weight(ifstep)                         
            wcvt(i,j)=wcvt(i,j)+weight(ifstep)                                  
          endif                                                                 
        enddo                                                                   
      enddo                                                                     
      do k=1,lsoil_                                                             
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            fsmc(i,j,k)=fsmc(i,j,k)+smc(i,j,k)*weight(ifstep)                   
            fstc(i,j,k)=fstc(i,j,k)+stc(i,j,k)*weight(ifstep)                   
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
      return                                                                    
      end                                                                       
