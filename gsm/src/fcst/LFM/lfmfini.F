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
      subroutine lfmfini(ipstep,fhour)                                          
c                                                                               
#include <paramodel.h>
#include <comfspec.vr.h>
#include <comfphys.h>
#include <comio.h>
#include <comlfm.h>
c                                                                               
      parameter(cvb0=100.,cvt0=0.)                                              
c                                                                               
c     dimension idate(4)                                                        
c                                                                               
#ifndef NOPRINT
      write(6,*) 'lfmfini called.  fhour=',fhour                                
#endif
c                                                                               
c  non-existent nlfmsgi is used for start of filtering integration              
c                                                                               
      read(nlfmsgi,end=999,err=999) lab                                         
      go to 888                                                                 
c                                                                               
c  start of new filtering integration                                           
c                                                                               
999   continue                                                                  
        ipstep=1                                                                
c                                                                               
c  zero out output surface array                                                
c                                                                               
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            islmsk(i,j,1)=0                                                     
            islmsk(i,j,2)=0                                                     
            islmsk(i,j,3)=0                                                     
          enddo                                                                 
        enddo                                                                   
        do l=1,lnt2_                                                            
          fq(l)=fq(l)+qm(l)*weight(1)                                           
        enddo                                                                   
        do k=1,levs_                                                            
          do l=1,lnt2_                                                          
            fte(l,k)=tem(l,k)*weight(1)                                         
            fdi(l,k)=dim(l,k)*weight(1)                                         
            fze(l,k)=zem(l,k)*weight(1)                                         
          enddo                                                                 
        enddo                                                                   
        do k=1,levh_                                                            
          do l=1,lnt2_                                                          
            frq(l,k)=rm(l,k)*weight(1)                                          
          enddo                                                                 
        enddo                                                                   
c                                                                               
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            ftsea(i,j)=tsea(i,j)*weight(1)                                      
            fsheleg(i,j)=sheleg(i,j)*weight(1)                                  
            ftg3(i,j)=tg3(i,j)*weight(1)                                        
            fzorl(i,j)=zorl(i,j)*weight(1)                                      
            fplantr(i,j)=plantr(i,j)*weight(1)                                  
            fcv(i,j)=cv(i,j)*weight(1)                                          
            do il = 1, 4
              falbedo(i,j,il)=albedo(i,j,il)*weight(1)                                  
            enddo
            ff10m(i,j)=f10m(i,j)*weight(1)                                      
            fcanopy(i,j)=fcanopy(i,j)*weight(1)                                 
            isl=nint(slmsk(i,j))+1                                              
            islmsk(i,j,isl)=islmsk(i,j,isl)+1                                   
            if(cvb(i,j).ne.cvb0) then                                           
              fcvb(i,j)=cvb(i,j)*weight(1)                                      
              wcvb(i,j)=weight(1)                                               
            else                                                                
              fcvb(i,j)=0.                                                      
              wcvb(i,j)=0.                                                      
            endif                                                               
            if(cvt(i,j).ne.cvt0) then                                           
              fcvt(i,j)=cvt(i,j)*weight(1)                                      
              wcvt(i,j)=weight(1)                                               
            else                                                                
              fcvt(i,j)=0.                                                      
              wcvt(i,j)=0.                                                      
            endif                                                               
          enddo                                                                 
        enddo                                                                   
        do k=1,lsoil_                                                           
          do j=1,latg2_                                                         
            do i=1,lonf2_                                                       
              fsmc(i,j,k)=smc(i,j,k)*weight(1)                                  
              fstc(i,j,k)=stc(i,j,k)*weight(1)                                  
            enddo                                                               
          enddo                                                                 
        enddo                                                                   
        return                                                                  
c                                                                               
c  continuation of filtering integration                                        
c                                                                               
  888 continue                                                                  
        read(nlfmsgi) ipstep,idate                                              
#ifndef NOPRINT
        write(6,*) 'ipstep,idate of filtered sig=',ipstep,idate                 
#endif
        read(nlfmsgi)                                                           
        read(nlfmsgi)(fq(i),i=1,lnt2_)                                          
        do k=1,levs_                                                            
          read(nlfmsgi) (fte(i,k),i=1,lnt2_)                                    
        enddo                                                                   
        do k=1,levs_                                                            
          read(nlfmsgi) (fdi(i,k),i=1,lnt2_)                                    
          read(nlfmsgi) (fze(i,k),i=1,lnt2_)                                    
        enddo                                                                   
        do k=1,levh_                                                            
          read(nlfmsgi) (frq(i,k),i=1,lnt2_)                                    
        enddo                                                                   
c                                                                               
        read(nlfmsfi) lab                                                       
        read(nlfmsfi) igstep,idate                                              
#ifndef NOPRINT
        write(6,*) 'igstep,idate of filtered sfc=',igstep,idate                 
#endif
        if(igstep.ne.ipstep) then                                               
          write(6,*) 'no. of steps on sig and sfc does not match'               
#ifdef MP
#ifdef RSM
          call rmpabort
#else
          call mpabort
#endif
#else
          call abort
#endif
        endif                                                                   
        read(nlfmsfi) ftsea                                                     
        read(nlfmsfi) fsmc                                                      
        read(nlfmsfi) fsheleg                                                   
        read(nlfmsfi) fstc                                                      
        read(nlfmsfi) ftg3                                                      
        read(nlfmsfi) fzorl                                                     
        read(nlfmsfi) fcv                                                       
        read(nlfmsfi) fcvb,wcvb                                                 
        read(nlfmsfi) fcvt,wcvt                                                 
        read(nlfmsfi) falbedo                                                   
        read(nlfmsfi) islmsk                                                    
        read(nlfmsfi) fplantr                                                   
        read(nlfmsfi) fcanopy                                                   
        read(nlfmsfi) ff10m                                                     
c                                                                               
      return                                                                    
      end                                                                       
