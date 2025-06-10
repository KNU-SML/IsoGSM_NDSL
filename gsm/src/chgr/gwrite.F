#define IGEN igen
#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
#include <define.h>
#ifndef HYBRID
      subroutine gwrite(n,fhour,idate,q,te,di,ze,chgrp,sl,si,gz,ltrn
#else
      subroutine gwrite(n,fhour,idate,q,te,di,ze,chgrp,ak5,bk5,gz,ltrn
#endif
     &       ,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
     &       ,iret)
#include <paramodel.h>
      save                                                                      
c                                                                               
#include <paramter.h>
c                                                                               
      parameter(idimt=idim*2,jdimhf=jdim/2,                                     
     1          mwavep=mwave+1,mwave2=mwave+2,                                  
     2          mdim=(mwave+1)*(mwave+2),                                       
     3          mdimv=(mwave+1)*(mwave+4),                                      
     4          mwvp2=(mwave+1)*2,                                              
     5          mdimhf=mdim/2,mdimvh=mdimv/2,                                   
     6          kdimp=kdim+1,kdimm=kdim-1)                                      
c                                                                               
#include <comcon.h>
c                                                                               
      common/comind/ indxnn(mdim),indxmm(mdim)                                  
      common /label/ lab(8)                                                     
      character*4 lab                                                           
c                                                                               
      dimension gz(mdim),q(mdim),                                               
     1          te(mdim,kdim),di(mdim,kdim),ze(mdim,kdim),                      
     2          chgrp(mdim,kdim)                                                
c                                                                               
#include <padchgr.h>
      parameter(npad=nscrch-ngwri)                                              
      common /scratch/ b(mdim),c(mdim),pad(npad)                                
c                                                                               
      dimension idate(4),si(kdimp),sl(kdim)                                     
      dimension ak5(kdimp),bk5(kdimp)
      dimension ak5x(kdimp),bk5x(kdimp)
c
#ifndef HYBRID
      parameter(kdum=201-levs_-1-levs_,kdum2=21,kens=2)
#else
      parameter(kdum=201-levs_-2-levs_,kdum2=21,kens=2)
#endif
      dimension dummy(kdum),dummy2(kdum2),ensemble(kens)
c                                                                               
      logical ltrn                                                              
c                                                                               
      call maklab(lab,idate,fhour,ifin)                                         
c                                                                               
      write(n) lab                                                              
c
c     write(n)fhour,idate,si,sl                                                 
c    &       ,dummy,waves,xlayers,trun,order,realform,gencode                   
c
      do i=1,kdum
        dummy(i)=0.
      enddo
      waves=jcap_
      xlayers=levs_
      trun=1.
      order=2.
      realform=1.
      gencode=igen_
      rlond=lonf_
      rlatd=latg_
      rlonp=lonf_
      rlatp=latg_
      rlonr=lonf_
      rlatr=latg_
      tracers=ntrac_
      clouds=ncldg_
      pdryini=0.
      subcen=0.
      do i=1,kens
        ensemble(i)=0.
      enddo
      ppid=0.
      slid=0.
      vcid=0.
      vmid=0.
      vtid=0.
      do k=1,kdum2
        dummy2(k)=0.
      enddo
c
#ifndef HYBRID
      write(n)fhour,idate,si,sl
#else
      do k=1,levs_+1
        ak5x(k)=ak5(k)*1000. ! cb -> Pa
        bk5x(k)=bk5(k)
      enddo
c
      write(n)fhour,idate,ak5x,bk5x
#endif
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
c                                                                               
      if(ltrn) then                                                             
      do 1010 i = 1, mdim                                                       
      b(indxmm(i)) = gz(i)                                                      
      c(indxmm(i)) =  q(i)                                                      
 1010 continue                                                                  
      write(n) b                                                                
      write(n) c                                                                
c                                                                               
      do 15 k=1,kdim                                                            
      do 1015 i = 1, mdim                                                       
      b(indxmm(i)) = te(i,k)                                                    
 1015 continue                                                                  
      write(n) b                                                                
   15 continue                                                                  
c                                                                               
      do 20 k=1,kdim                                                            
      do 1020 i = 1, mdim                                                       
      b(indxmm(i)) = di(i,k)                                                    
      c(indxmm(i)) = ze(i,k)                                                    
 1020 continue                                                                  
      write(n) b                                                                
      write(n) c                                                                
   20 continue                                                                  
c                                                                               
      do 25 k=1,kdimq                                                           
      do 1025 i = 1, mdim                                                       
      b(indxmm(i)) = chgrp(i,k)                                                 
 1025 continue                                                                  
      write(n) b                                                                
   25 continue                                                                  
c                                                                               
      else                                                                      
c                                                                               
      write(n) gz                                                               
      write(n) q                                                                
c                                                                               
      do 55 k=1,kdim                                                            
      write(n) (te(i,k),i=1,mdim)                                               
   55 continue                                                                  
c                                                                               
      do 60 k=1,kdim                                                            
      write(n) (di(i,k),i=1,mdim)                                               
      write(n) (ze(i,k),i=1,mdim)                                               
   60 continue                                                                  
c                                                                               
      do 65 k=1,kdimq                                                           
      write(n) (chgrp(i,k),i=1,mdim)                                            
   65 continue                                                                  
c                                                                               
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
