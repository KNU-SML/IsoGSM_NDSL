#define CAL cal
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LTSTP ltstp
#define MLVARK mlvark
#define NCLDG ncldg
#define NTRAC ntrac
#define SLVARK slvark

      subroutine kenput(lat,rcl,                                                
     &                  slmsk,psexp,tg3,sheleg,radsl,dlwsf1,                    
     &                  tsea,qss,plantr,gflx,zorl,cd,cdq,                       
     &                  rnet,hflx,stsoil,                                  
     &                  canopy,drain,smsoil,runof,cld1d,                        
     &                  u10,v10,t2,q2,                                          
     &                  hpbl,gamt,gamq,                                         
     &                  dqsfc1,dtsfc1,dusfc1,dvsfc1,                            
     &                  dusfcg,dvsfcg,                                          
     &                  rainc,rainl,                                            
     &                  u,v,t,q,ntotal,hsw,hlw,vvel,qcicnv,qrscnv,
     &                  snowmt,snowev,snowfl)                                   
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
#include <comgpd.h>
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
      parameter(im=LONF2S,ix=LONF2S,km=levs_) 
      real slmsk(im),psexp(im),tg3(im),sheleg(im),radsl(im),dlwsf1(im),         
     &     tsea(im),qss(im),plantr(im),gflx(im),zorl(im),cd(im),cdq(im),        
     &     rnet(im),hflx(im),stsoil(im,2),                             
     &     canopy(im),drain(im),smsoil(im,2),runof(im),cld1d(im),               
     &     u10(im),v10(im),t2(im),q2(im),                                       
     &     hpbl(im),gamt(im),gamq(im),                                          
     &     dqsfc1(im),dtsfc1(im),dusfc1(im),dvsfc1(im),                         
     &     dusfcg(im),dvsfcg(im),                                               
     &     rainc(im),rainl(im),                                                 
     &     u(ix,km),v(ix,km),t(ix,km),q(ix,km,ntotal),                                 
     &     hsw(im,km),hlw(im,km),vvel(im,km),
     &     qcicnv(im,km),qrscnv(im,km),                                   
     &     snowmt(im),snowev(im),snowfl(im)
      parameter(cnwatt=-cal_*1.e4/60.)                                         
#ifdef KEN
      do 330 igpt=1,npoint                                                      
      if(lat.eq.jgrd(igpt)) then                                                
      svdata(11,igpt,itnum)=svdata(11,igpt,itnum)+rainc(igrd(igpt))             
      svdata(12,igpt,itnum)=svdata(12,igpt,itnum)+rainl(igrd(igpt))             
      if(isave.ne.0) then                                                       
      svdata(  1,igpt,itnum)= igrd(igpt)                                        
      svdata(  2,igpt,itnum)= jgrd(igpt)                                        
      svdata(  3,igpt,itnum)= slmsk (igrd(igpt))                                
      svdata(  4,igpt,itnum)= psexp (igrd(igpt)) *10.                           
      svdata(  8,igpt,itnum)= tg3 (igrd(igpt))                                  
      svdata( 10,igpt,itnum)= sheleg(igrd(igpt))                                
      svdata( 13,igpt,itnum)= radsl(igrd(igpt))*cnwatt                          
      svdata( 14,igpt,itnum)= dlwsf1(igrd(igpt))                                
      svdata(  5,igpt,itnum)= tsea (igrd(igpt))                                 
      svdata( 15,igpt,itnum)= qss  (igrd(igpt))                                 
      svdata( 16,igpt,itnum)= plantr(igrd(igpt))                                
      svdata( 19,igpt,itnum)= gflx(igrd(igpt))                                  
      svdata( 22,igpt,itnum)= zorl (igrd(igpt))                                 
      svdata( 23,igpt,itnum)= cd  (igrd(igpt))                                  
      svdata( 24,igpt,itnum)= cdq (igrd(igpt))                                  
      svdata( 62,igpt,itnum)= rnet (igrd(igpt))                                 
#endif
ccken svdata( 63,igpt,itnum)= evap (igrd(igpt))                                 
#ifdef KEN
      svdata( 64,igpt,itnum)= hflx (igrd(igpt))                                 
      svdata(  6,igpt,itnum)= stsoil (igrd(igpt),1)                             
      svdata(  7,igpt,itnum)= stsoil (igrd(igpt),2)                             
      svdata( 34,igpt,itnum)= u10 (igrd(igpt))                                  
      svdata( 35,igpt,itnum)= v10 (igrd(igpt))                                  
      svdata( 30,igpt,itnum)= t2  (igrd(igpt))                                  
      svdata( 31,igpt,itnum)= q2  (igrd(igpt))                                  
      svdata( 32,igpt,itnum)= canopy(igrd(igpt))                                
      svdata( 33,igpt,itnum)= drain(igrd(igpt))                                 
      svdata( 17,igpt,itnum)= dqsfc1(igrd(igpt))                                
      svdata( 18,igpt,itnum)= dtsfc1(igrd(igpt))                                
      svdata( 20,igpt,itnum)= dusfc1(igrd(igpt))                                
      svdata( 21,igpt,itnum)= dvsfc1(igrd(igpt))                                
      svdata( 28,igpt,itnum)= dusfcg(igrd(igpt))                                
      svdata( 29,igpt,itnum)= dvsfcg(igrd(igpt))                                
      svdata(  9,igpt,itnum)= smsoil (igrd(igpt),1)                             
      svdata( 61,igpt,itnum)= smsoil (igrd(igpt),2)                             
#endif
ccken svdata( 65,igpt,itnum)= runof (igrd(igpt))                                
#ifdef KEN
      svdata( 63,igpt,itnum)= hpbl(igrd(igpt))                                  
      svdata( 64,igpt,itnum)= gamt(igrd(igpt))                                  
      svdata( 65,igpt,itnum)= gamq(igrd(igpt))                                  
      svdata( 66,igpt,itnum)= cld1d(igrd(igpt))                                
      svdata( 67,igpt,itnum)= snowev(igrd(igpt))                                
      svdata( 68,igpt,itnum)= snowmt(igrd(igpt))                                
      svdata( 69,igpt,itnum)= snowfl(igrd(igpt))                                
      if(ilshrt.lt.2) then                                                      
      r=sqrt(rcl)                                                               
      do 331 k=1,levs_                                                          
      svdata(k+slvark_+0*levs_,igpt,itnum)=u(igrd(igpt),k)*r                    
      svdata(k+slvark_+1*levs_,igpt,itnum)=v(igrd(igpt),k)*r                    
      svdata(k+slvark_+2*levs_,igpt,itnum)=t(igrd(igpt),k)                      
      svdata(k+slvark_+3*levs_,igpt,itnum)=q(igrd(igpt),k)                      
      if(ilshrt.lt.1) then                                                      
      svdata(k+slvark_+4*levs_,igpt,itnum)=hsw(igrd(igpt),k)                    
      svdata(k+slvark_+5*levs_,igpt,itnum)=hlw(igrd(igpt),k)                    
      svdata(k+slvark_+6*levs_,igpt,itnum)=vvel(igrd(igpt),k)                   
      if(ntrac_.ge.2) then
      svdata(k+slvark_+7*levs_,igpt,itnum)=q(igrd(igpt),k,2)
      endif
      svdata(k+slvark_+9*levs_,igpt,itnum)=qcicnv(igrd(igpt),k)
      svdata(k+slvark_+10*levs_,igpt,itnum)=qrscnv(igrd(igpt),k)
      if(ncldg_.ge.1) then
      do icloud = 1,ncldg_
      ic = ntrac_ + icloud
      ivar = 10 + icloud
      svdata(k+slvark_+ivar*levs_,igpt,itnum)=q(igrd(igpt),k,ic)
      enddo
      endif
      endif                                                                     
#endif
331   continue                                                                  
#ifdef KEN
      endif                                                                     
      endif                                                                     
      endif                                                                     
#endif
330   continue                                                                  
      return                                                                    
      end                                                                       
