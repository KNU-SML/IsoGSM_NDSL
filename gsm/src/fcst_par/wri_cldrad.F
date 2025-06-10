#define JCAP jcap
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LNT2 lnt2
#define LSOIL lsoil
#define NPES npes
#include <define.h>
      subroutine wri_cldrad(fhour,thour,idate,sl,colrad,
     2     fluxr,cvavg,qcicnv,qrscnv, qcilrg,qrslrg,taucld,cldwp,cldip,
     3     ndg)                                
c======================================================================
c     to write diagnostic data concerning with cloud & radiation
c     y.-h. byun               18 Aug 2004
c======================================================================
#include <paramodel.h>
c                                                                               
#ifdef CRAY_THREAD
cfpp$ expand(znlacm,accdia)
#endif
#ifndef RSM
#ifdef MP
#define LONF2S lonf2p_
#define LONF22S lonf22p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LONF22S lonf22_
#define LATG2S latg2_
#endif
#else
#ifdef MP
#ifdef RSMVECTORIZE
#define LONF2S igrd12p_*jgrd12p_
#else
#define LONF2S igrd12p_
#endif
#define LATG2S jgrd12p_
#else
#undef RSMVECTORIZE
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif
#ifndef RSM		/* not RSM */
#ifdef MP
#define MPGP2F mpgp2f
#define MPABORT mpabort
#else                   /* not MP */
#define MPABORT abort
#endif                  /* MP */
#else                   /* RSM */
#ifdef MP		/* MP */
#define MPGP2F rmpgp2f
#define MPABORT rmpabort
#include <comrmpi.h>
#else                   /* not MP */
#define MPABORT abort
#endif                  /* MP */
#endif                  /* RSM */
#ifdef DG3
      integer idate(4)                                                          
      dimension sl(levs_),colrad(latg2_)                                        
#endif
#ifdef MP
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
cld force                                                                       
#ifdef CLR
      dimension cfsw(lonf2_,latg2_,3)
      dimension cflw(lonf2_,latg2_,3)
      equivalence (cfsw(1,1,1),cflw(1,1,1))                                     
#endif
cld force                                                                       
#ifdef DG3
      dimension   fluxr(LONF2S,LATG2S,26)                 
      dimension   cvavg(LONF2S,LATG2S)                   
#endif
#ifdef CONDENS
      dimension
     + qcicnv(LONF2S,levs_,LATG2S), qrscnv(LONF2S,levs_,LATG2S),
     + qcilrg(LONF2S,levs_,LATG2S), qrslrg(LONF2S,levs_,LATG2S),
     + taucld(LONF2S,levs_,LATG2S),
     + cldwp(LONF2S,levs_,LATG2S),  cldip(LONF2S,levs_,LATG2S)
#endif
c                                                                               
#include <comio.h>
c                                                                               
#include <comgda.h>
#ifdef DG3
      dimension gda(nwgda)                                                      
#endif
c                                                                               
      parameter(iprs=1,itemp=11,iznlw=33,imerw=34,isphum=51,                    
     $          ipcpr=59,isnowd=65,icldf=71,iccldf=72,                          
     $          islmsk=81,izorl=83,ialbdo=84,isoilm=144,icemsk=91,              
     $          ilhflx=121,ishflx=122,izws=124,imws=125,ighflx=155,             
     $          iuswfc=160,idswfc=161,iulwfc=162,idlwfc=163,                    
     $          inswfc=164,inlwfc=165,                                          
     $          idswvb=166,idswvd=167,idswnb=168,idswnd=169,                    
     $          isglyr=175,icnpy=145,                                           
     $          idswf=204,idlwf=205,iuswf=211,iulwf=212,icpcpr=214)             
      parameter(isfc=1,itoa=8,ielev=105,                                        
     $          isglev=107,idbls=111,i2dbls=112,icolmn=200,                     
     $          ilcbl=212,ilctl=213,ilclyr=214,                                 
     $          imcbl=222,imctl=223,imclyr=224,                                 
     $          ihcbl=232,ihctl=233,ihclyr=234)                                 
      parameter(inst=10,iavg=3,iacc=4)                                          
      parameter(ifhour=1,ifday=2)                                               
      parameter(lonb=lonf2_/2,latb=latg2_*2)
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#else
#include <commpi.h>
#endif
#endif
#ifdef DG3
      logical lbm(lonf2_,latg2_)                                                
      character g(200+lonf2_*latg2_*(32+1)/8)                                   
      dimension ids(255)                                                        
      dimension iens(5)                                                         
      dimension iclyr(3),ictl(3),icbl(3),itlcf(3)                               
      data iclyr/ihclyr,imclyr,ilclyr/                                          
      data ictl /ihctl ,imctl ,ilctl /                                          
      data icbl /ihcbl ,imcbl ,ilcbl /                                          
      data itlcf/itoa,isfc,icolmn/                                              
      save phour
      data phour/0.0/
c
#endif
#ifdef ASSIGN
      character*120 asgnstr
#endif
c                                                                               
#ifdef DG3
      dimension work(lonf2_,latg2_),slmsep(lonf2_,latg2_),                      
     1          work2(lonf2_,latg2_,levs_+2/levs_),                             
     1          workc(lonf2_,20),iworkc(lonf2_)                                 
      dimension fluxw(lonf2_,latg2_,26)
#ifdef MP
      dimension work2p(lonf2p_,latg2p_,levs_+2/levs_)
#endif
c
      character*80 fni
      parameter(fni='clddiag',nchi=7)
      character*80 fno
c                                                                               
#ifdef SCM
      cl1=90.0
      idrt=0
#else
      cl1=colrad(1)
      idrt=4
#endif
c
#ifdef MP
      if(iope) then
        print *,' start wridia '
#endif
        call fnam(fni,nchi,thour,fno,ncho)
#ifdef ASSIGN
        write(asgnstr,'(23hassign -s unblocked  u:,I2,)') ndg
        call assign('assign -R')
        call assign(asgnstr)
#endif
        open(unit=ndg,file=fno(1:ncho),form='unformatted',err=900)
        go to 901
  900   continue
        write(6,*) ' error in opening file ',fno(1:ncho)
	call MPABORT
  901   continue
#ifndef NOPRINT
        write(6,*) ' file ',fno(1:ncho),' opened. unit=',ndg
#endif
        rewind ndg
#ifdef MP
      endif	! iope
#endif
      if( phour.eq.0.0 ) phour=fhour
      call idsdef(1,ids)                                                        
      iens(1)=1                                                                 
      iens(2)=ienst                                                             
      iens(3)=iensi                                                             
      iens(4)=1                                                                 
      iens(5)=255                                                               
      iyr=idate(4)                                                              
      imo=idate(2)                                                              
      ida=idate(3)                                                              
      ihr=idate(1)                                                              
      iftime=ifhour                                                             
      ifhr=nint(phour)                                                          
      ithr=nint(thour)                                                          
      dhour=thour-phour                                                         
      if(dhour.gt.0.) then                                                      
        rtime=1./(3600.*dhour)                                                  
      else                                                                      
        rtime=0.                                                                
      endif                                                                     
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(dhour.gt.0.) then                                                      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
#ifdef MP
      call MPGP2F(fluxr,lonf2p_,latg2p_,fluxw,lonf2_,latg2_,26)
      if(iope) then
#else
      do k=1,26
       do j=1,latg2_
        do i=1,lonf2_
          fluxw(i,j,k)=fluxr(i,j,k)
        enddo
       enddo
      enddo
#endif
#ifdef MP
      endif ! iope
      call MPGP2F(cvavg,lonf2p_,latg2p_,work,lonf2_,latg2_,1)
      if(iope) then
#endif
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
#ifdef MP
          work(i,j)= work(i,j)*rtime*1.e2                                             
#else
          work(i,j)=cvavg(i,j)*rtime*1.e2                                       
#endif
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,iccldf,isfc,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iccldf),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c
#ifdef CLR
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          work(i,j)=fluxw(i,j,21)
          work(i,j)=work(i,j)*rtime                                             
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,iulwfc,itoa,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iulwfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          work(i,j)=fluxw(i,j,22)
          work(i,j)=work(i,j)*rtime                                             
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,iuswfc,itoa,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          work(i,j)=fluxw(i,j,25)
          work(i,j)=work(i,j)*rtime                                             
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,idlwfc,isfc,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idlwfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          work(i,j)=fluxw(i,j,23)
          work(i,j)=work(i,j)*rtime                                             
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,idswfc,isfc,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c
      do j=1,latg2_                                                             
        do i=1,lonf2_                                                           
          work(i,j)=fluxw(i,j,24)
          work(i,j)=work(i,j)*rtime                                             
        enddo                                                                   
      enddo                                                                     
      call rowsep(work)                                                         
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,iuswfc,isfc,0,0,iyr,imo,ida,ihr,                            
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
c.....   compute sw cloud forcing at toa (cld-clear),                           
c         flip sign so positive means cld is warming relative to clear          
      do 600 j=1,latg2_                                                         
       do 600 i=1,lonf2_                                                        
        cfsw(i,j,1) = - (fluxw(i,j,2)-fluxw(i,j,22))                            
  600 continue                                                                  
c.....    compute cloud forcing at sfc (cld-clear)                              
c    again flip sign so positive means cld is warming relative to clear         
      do 601 j=1,latg2_                                                         
       do 601 i=1,lonf2_                                                        
       cfsw(i,j,2) =-(fluxw(i,j,3)-fluxw(i,j,4)                                 
     1                      -(fluxw(i,j,24)-fluxw(i,j,23)))                     
  601 continue                                                                  
c.....  flip sign so positive means cld is warming relative to clear            
      do 602 j=1,latg2_                                                         
       do 602 i=1,lonf2_                                                        
        cfsw(i,j,3) = - (cfsw(i,j,2) - cfsw(i,j,1))                             
  602 continue                                                                  
      do 603 k=1,3                                                              
c
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            work(i,j)=cfsw(i,j,k)
            work(i,j)=work(i,j)*rtime                                           
          enddo                                                                 
        enddo                                                                   
       call rowsep(work)                                                        
cc     lbm=work.ne.0.                                                           
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,inswfc,itlcf(k),0,0,iyr,imo,ida,ihr,                        
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(inswfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
  603 continue                                                                  
c.....   compute lw cloud forcing at toa (cld-clear),                           
c         flip sign so positive means cld is warming relative to clear          
      do 1600 j=1,latg2_                                                        
       do 1600 i=1,lonf2_                                                       
        cflw(i,j,1) = - (fluxw(i,j,1)-fluxw(i,j,21))                            
 1600 continue                                                                  
c.....    compute cloud forcing at sfc (cld-clear)                              
c    again flip sign so positive means cld is warming relative to clear         
      do 1601 j=1,latg2_                                                        
       do 1601 i=1,lonf2_                                                       
       cflw(i,j,2) = - (fluxw(i,j,25)-fluxw(i,j,19))                            
 1601 continue                                                                  
c.....  flip sign so positive means cld is warming relative to clear            
      do 1602 j=1,latg2_                                                        
       do 1602 i=1,lonf2_                                                       
        cflw(i,j,3) = - (cflw(i,j,2) - cflw(i,j,1))                             
 1602 continue                                                                  
      do 1603 k=1,3                                                             
c
       do j=1,latg2_                                                            
         do i=1,lonf2_                                                          
           work(i,j)=cflw(i,j,k)
           work(i,j)=work(i,j)*rtime                                            
         enddo                                                                  
       enddo                                                                    
       call rowsep(work)                                                        
      call gribit(work,lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,                
     &            0,inlwfc,itlcf(k),0,0,iyr,imo,ida,ihr,                        
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(inlwfc),iens,             
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
 1603 continue                                                                  
#endif
cld force                                                                       
c..                                                                             
csib                                                                            
ccdg3 work=work*rtime                                                           
#ifdef MP
      endif 	! iope
#endif
c
      do 750 kd=1,kdgda                                                         
c
      do 720 j=1,nrgda                                                          
      kgda=(j-1)*kdgda+kd                                                       
      call getdia(kgda,nwgda,gda)                                               
      igda=0                                                                    
      do 710 k=1,levs_                                                          
      do 710 i=1,LONF2S                                                         
      igda=igda+1                                                               
#ifdef MP
      work2p(i,j,k)=gda(igda)*rtime
#else
      work2(i,j,k)=gda(igda)*rtime                                              
#endif
710   continue                                                                  
720   continue                                                                  
#ifdef MP
      call MPGP2F(work2p,lonf2p_,latg2p_,work2,lonf2_,latg2_,levs_)
      if(iope) then
#endif
      ipu=ipugda(kd)                                                            
      ibm=ibmgda(kd)                                                            
      do 730 k=1,levs_                                                          
      isl=nint(sl(k)*1.e4)                                                      
      call rowsep(work2(1,1,k))                                                 
      ibitmap=0
      if(ibm.ne.0) then                                                         
        ibitmap=1
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            lbm(i,j)=work2(i,j,k).ne.0.                                         
          enddo                                                                 
        enddo                                                                   
      endif                                                                     
      call gribit(work2(1,1,k),lbm,idrt,lonb,latb,16,cl1,28,2,icen,igen,        
c    &            ibitmap,ipu,isglev,0,isl,iyr,imo,ida,ihr,
     &                  0,ipu,isglev,0,isl,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(ipu),iens,                
     &            0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)                            
      if(ierr.eq.0) call wryte(ndg,lg,g)                                        
730   continue                                                                  
c                                                                               
#ifdef MP
      endif ! iope
#endif
750   continue                                                                  
cyhb
c
c ... write more ...
c
      do nvar = 1, 7
       do k=1,levs_
       do j=1,LATG2S
       do i=1,LONF2S
#ifdef MP
         if (nvar.eq.1) work2p(i,j,k)=qcicnv(i,k,j)
         if (nvar.eq.2) work2p(i,j,k)=qrscnv(i,k,j)
         if (nvar.eq.3) work2p(i,j,k)=qcilrg(i,k,j)
         if (nvar.eq.4) work2p(i,j,k)=qrslrg(i,k,j)
         if (nvar.eq.5) work2p(i,j,k)=taucld(i,k,j)
         if (nvar.eq.6) work2p(i,j,k)=cldwp (i,k,j)
         if (nvar.eq.7) work2p(i,j,k)=cldip (i,k,j)
#else
         if (nvar.eq.1) work2(i,j,k)=qcicnv(i,k,j)
         if (nvar.eq.2) work2(i,j,k)=qrscnv(i,k,j)
         if (nvar.eq.3) work2(i,j,k)=qcilrg(i,k,j)
         if (nvar.eq.4) work2(i,j,k)=qrslrg(i,k,j)
         if (nvar.eq.5) work2(i,j,k)=taucld(i,k,j)
         if (nvar.eq.6) work2(i,j,k)=cldwp (i,k,j)
         if (nvar.eq.7) work2(i,j,k)=cldip (i,k,j)
#endif
       enddo
       enddo
       enddo
#ifdef MP
       call MPGP2F(work2p,lonf2p_,latg2p_,work2,lonf2_,latg2_,levs_)
       if(iope) then
#endif
       ipunv = 185 + nvar - 1
       idsnv = 3
       if (nvar.le.4) then
         idsnv = 4
       elseif (nvar.ge.6) then
         idsnv = 1
       endif
       do k=1,levs_
        isl=nint(sl(k)*1.e4)
        call rowsep(work2(1,1,k))
        call gribit(work2(1,1,k),lbm,idrt,lonb,latb,16,cl1,28,2,
     +              icen,igen,
     &              0,ipunv,isglev,0,isl,iyr,imo,ida,ihr,
     &              iftime,ifhr,ithr,iavg,0,0,icen2,idsnv,iens,
     &              0.,0.,0.,0.,0.,0.,0.,0.,g,lg,ierr)
        if(ierr.eq.0) call wryte(ndg,lg,g)
       enddo
#ifdef MP
       endif ! iope
#endif
      enddo
cyhb
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif                                                                     
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c                                                                               
      phour=thour
c                                                                               
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(ndg)                                                                
#ifdef MP
      endif
#endif
c 
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#endif	
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c                                                                               
      return                                                                    
      end                                                                       
