#define CP cp
#define G g
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define OMEGA omega
#define RD rd
#define RERTH rerth
#define RV rv
#include <define.h>
      subroutine znlavb(lat,lonb2,lonb22,levs,secphy,secrad,                    
     &                  rbs2,si,sl,del,psexp,geshem,bengsh,                     
     &                  dusfc,dvsfc,dtsfc,dqsfc,dugwd,dvgwd,                    
     &                  slmsk,sheleg,tsea,soilm,tg1,tg2,tg3,                    
     &                  zorl,ep,cldwrk,dlwsfc,ulwsfc,                           
     &                  dg,tg,rqg,dlam,dphi,ug,vg)                              
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
cfpp$ expand(fpvs,fpvs0,znlaim,znlais)                                                
#endif
c                                                                               
#include <paramodel.h>
      dimension si(lonb2,levs+1),sl(lonb2,levs),del(lonb2,levs)                                   
      dimension psexp(lonb2)                                                    
      dimension geshem(lonb2),bengsh(lonb2)                                     
      dimension dusfc(lonb2),dvsfc(lonb2)                                       
      dimension dtsfc(lonb2),dqsfc(lonb2)                                       
      dimension dugwd(lonb2),dvgwd(lonb2)                                       
      dimension slmsk(lonb2),sheleg(lonb2)                                      
      dimension tsea(lonb2),soilm(lonb2)                                        
      dimension tg1(lonb2),tg2(lonb2),tg3(lonb2)                                
      dimension zorl(lonb2),ep(lonb2),cldwrk(lonb2)                             
      dimension dlwsfc(lonb2),ulwsfc(lonb2)                                     
      dimension dg(lonb22,levs)                                                 
      dimension tg(lonb22,levs),rqg(lonb22,levs)                                
      dimension dlam(lonb2),dphi(lonb2)                                         
      dimension ug(lonb22,levs),vg(lonb22,levs)                                 
c                                                                               
#include <comznl.h>
c                                                                               
      parameter(cp=cp_,rd=rd_,rv=rv_)                                           
      parameter(fkappa=rd/cp,fvirt=rv/rd-1.,feps=rd/rv,fepsm1=rd/rv-1.)         
#ifdef DYNAMIC_ALLOC
      dimension ist(lonb2)                                                      
      dimension work(lonb2,levs),work1(lonb2)                                   
#else
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
      dimension ist(LONF2S)
      dimension work(LONF2S,levs_),work1(LONF2S)
#endif
                                                                                
c                                                                               
      lonb=lonb2/2                                                              
c                                                                               
      zwm(1,lat)=0.                                                             
      zwm(2,lat)=0.                                                             
      do i=1,lonb                                                               
        zwm(1,lat)=zwm(1,lat)+psexp(i)                                          
        zwm(2,lat)=zwm(2,lat)+psexp(i+lonb)                                     
      enddo                                                                     
c                                                                               
      do i=1,lonb2                                                              
        work1(i)=0.                                                             
      enddo                                                                     
      do k=levs,1,-1                                                            
        do i=1,lonb2                                                            
          cg1=rbs2*(ug(i,k)*dlam(i)+vg(i,k)*dphi(i))                            
          cg2=cg1+dg(i,k)                                                       
          work1(i)=work1(i)+cg2*(sl(i,k)-si(i,k+1))                                 
          work(i,k)=fkappa*tg(i,k)*(cg1-work1(i)/sl(i,k))                         
          work1(i)=work1(i)+cg2*(si(i,k)-sl(i,k))                                   
        enddo                                                                   
      enddo                                                                     
      call znlaim(lonb2,work,psexp,nmomega,lat)                                 
c                                                                               
      do k=1,levs                                                               
        do i=1,lonb2                                                            
          work(i,k)=tg(i,k)*(1.+fvirt*rqg(i,k))                                 
        enddo                                                                   
      enddo                                                                     
      call znlaim(lonb2,work,psexp,nmtv,lat)                                    
c                                                                               
      do k=1,levs                                                               
        do i=1,lonb2                                                            
#ifdef ICE
          es=fpvs(tg(i,k))                                                      
#else
          es=fpvs0(tg(i,k))                                                      
#endif
          qs=feps*es/(sl(i,k)*psexp(i)+fepsm1*es)                                 
          work(i,k)=100.*rqg(i,k)/qs                                            
        enddo                                                                   
      enddo                                                                     
      call znlaim(lonb2,work,psexp,nmrh,lat)                                    
c                                                                               
      do k=1,levs                                                               
        do i=1,lonb2                                                            
          work(i,k)=rbs2*0.5*(ug(i,k)**2+vg(i,k)**2)                            
        enddo                                                                   
      enddo                                                                     
      call znlaim(lonb2,work,psexp,nmke,lat)                                    
c                                                                               
      do k=1,levs                                                               
        do i=1,lonb2                                                            
          work(i,k)=dg(i,k)**2                                                  
        enddo                                                                   
      enddo                                                                     
      call znlaim(lonb2,work,psexp,nmdiv2,lat)                                  
c                                                                               
      call znlaim(lonb22,ug,psexp,nmu,lat)                                      
      call znlaim(lonb22,vg,psexp,nmv,lat)                                      
      call znlaim(lonb22,tg,psexp,nmt,lat)                                      
      call znlaim(lonb22,rqg,psexp,nmq,lat)                                     
c                                                                               
      wginst=1.                                                                 
      wgtime=0.                                                                 
      if(secphy.gt.0.) wgtime=1./secphy                                         
      wrtime=0.                                                                 
      if(secrad.gt.0.) wrtime=1./secrad                                         
      zhm(nmu)=wginst                                                           
      zhm(nmv)=wginst                                                           
      zhm(nmtv)=wginst                                                          
      zhm(nmq)=wginst                                                           
      zhm(nmvot2)=0.                                                            
      zhm(nmdiv2)=wginst                                                        
      zhm(nmomega)=wginst                                                       
      zhm(nmt)=wginst                                                           
      zhm(nmrh)=wginst                                                          
      zhm(nmke)=wginst                                                          
      zhm(nmtconv)=wgtime                                                       
      zhm(nmtlarg)=wgtime                                                       
      zhm(nmtshal)=wgtime                                                       
      zhm(nmtvrdf)=wgtime                                                       
      zhm(nmqconv)=wgtime                                                       
      zhm(nmqshal)=wgtime                                                       
      zhm(nmqvrdf)=wgtime                                                       
      zhm(nmuvrdf)=wgtime                                                       
      zhm(nmvvrdf)=wgtime                                                       
      zhm(nmthsw)=wgtime                                                        
      zhm(nmthlw)=wgtime                                                        
      zhm(nmtcld)=wrtime                                                        
      zhm(nmtccv)=wrtime                                                        
c                                                                               
      do i=1,lonb2                                                              
        if(slmsk(i).eq.0.) then                                                 
          ist(i)=6                                                              
        elseif(slmsk(i).eq.1.) then                                             
          ist(i)=2                                                              
          if(sheleg(i).gt.1.e-3) ist(i)=3                                       
        else                                                                    
          ist(i)=4                                                              
          if(sheleg(i).gt.1.e-3) ist(i)=5                                       
        endif                                                                   
      enddo                                                                     
      zws(1,1,lat)=lonb                                                         
      zws(2,1,lat)=lonb                                                         
      do k=2,nst                                                                
        zws(1,k,lat)=0.                                                         
        zws(2,k,lat)=0.                                                         
      enddo                                                                     
      do i=1,lonb                                                               
        zws(1,ist(i),lat)=zws(1,ist(i),lat)+1.                                  
        zws(2,ist(i+lonb),lat)=zws(2,ist(i+lonb),lat)+1.                        
      enddo                                                                     
c                                                                               
      acl=rerth_/sqrt(rbs2)                                                     
      do i=1,lonb2                                                              
        work1(i)=0.                                                             
      enddo                                                                     
      do k=1,levs                                                               
        do i=1,lonb2                                                            
          work1(i)=work1(i)+del(i,k)*ug(i,k)                                      
        enddo                                                                   
      enddo                                                                     
      uafac=rerth_*1.e3/g_                                                      
      do i=1,lonb2                                                              
        work1(i)=uafac*psexp(i)*work1(i)                                        
      enddo                                                                     
      call znlais(work1,ist,1.,nsua,lat)                                        
c                                                                               
      uafac=omega_*acl*acl*1.e3/g_                                              
      uaref=omega_*acl*acl*1.e4                                                 
      do i=1,lonb2                                                              
        work1(i)=uafac*psexp(i)-uaref                                           
      enddo                                                                     
      call znlais(work1,ist,1.,nsuap,lat)                                       
c                                                                               
      do i=1,lonb2                                                              
        work1(i)=tg(i,1)*(1.+fvirt*rqg(i,1))                                    
      enddo                                                                     
      call znlais(work1,ist,1.,nstvs,lat)                                       
c                                                                               
      do i=1,lonb2                                                              
#ifdef ICE
        es=fpvs(tg(i,1))                                                        
#else
        es=fpvs0(tg(i,1))                                                        
#endif
        qs=feps*es/(sl(i,1)*psexp(i)+fepsm1*es)                                   
        work1(i)=100.*rqg(i,1)/qs                                               
      enddo                                                                     
      call znlais(work1,ist,1.,nsrhs,lat)                                       
c                                                                               
      call znlais(sheleg,ist,1.,nssnow,lat)                                     
      call znlais(psexp,ist,1.,nsps,lat)                                        
      call znlais(tg(1,1),ist,1.,nsts,lat)                                      
      call znlais(rqg(1,1),ist,1.,nsqs,lat)                                     
      call znlais(tsea,ist,1.,nstskin,lat)                                      
      call znlais(soilm,ist,1.,nswet,lat)                                       
      call znlais(tg1,ist,1.,nstg1,lat)                                         
      call znlais(tg2,ist,1.,nstg2,lat)                                         
      call znlais(tg3,ist,1.,nstg3,lat)                                         
      call znlais(zorl,ist,1.,nszorl,lat)                                       
c                                                                               
      if(secphy.gt.0.) then                                                     
        do i=1,lonb2                                                            
          work1(i)=0.                                                           
          if(geshem(i).gt.0.) work1(i)=100.                                     
        enddo                                                                   
        call znlais(work1,ist,1.,nsrcov,lat)                                    
c                                                                               
        do i=1,lonb2                                                            
          work1(i)=0.                                                           
          if(bengsh(i).gt.0.) work1(i)=100.                                     
        enddo                                                                   
        call znlais(work1,ist,1.,nsrcovc,lat)                                   
c                                                                               
        do i=1,lonb2                                                            
          work1(i)=ulwsfc(i)-dlwsfc(i)                                          
        enddo                                                                   
        call znlais(work1,ist,1./secphy,nssfclw,lat)                            
c                                                                               
        call znlais(geshem,ist,86400./secphy,nsrain,lat)                        
        call znlais(bengsh,ist,86400./secphy,nsrainc,lat)                       
        call znlais(dusfc,ist,1./secphy,nsusfc,lat)                             
        call znlais(dvsfc,ist,1./secphy,nsvsfc,lat)                             
        call znlais(dugwd,ist,1./secphy,nsugwd,lat)                             
        call znlais(dvgwd,ist,1./secphy,nsvgwd,lat)                             
        call znlais(dusfc,ist,acl/secphy,nsuasfc,lat)                           
        call znlais(dugwd,ist,acl/secphy,nsuagwd,lat)                           
        call znlais(dtsfc,ist,1./secphy,nstsfc,lat)                             
        call znlais(dqsfc,ist,1./secphy,nsqsfc,lat)                             
        call znlais(ep,ist,1./secphy,nsep,lat)                                  
        call znlais(cldwrk,ist,1./secphy,nscldwrk,lat)                          
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
