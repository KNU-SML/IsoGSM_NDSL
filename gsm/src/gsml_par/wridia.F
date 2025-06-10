#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNTP lntp
#define LNUV lnuv
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NTOTAL ntotal
#define NTRAC ntrac
      subroutine wridia(slmask,ndg)
c
#include <paramodel.h>
#ifndef RSM
#include <comfphys.h>
#include <comfspec.vr.h>
#include <comfver.h>
#include <comfgrid.h>
#include <radiag.h>
#define LONF2 lonf2_
#define LATG2 latg2_
#define ILEN lonf_
#define JLEN latg_
#ifdef MP
#define MPGP2F mpgp2f
#define MPABORT mpabort
#include <commpi.h>
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else			/* MP */
#define MPABORT abort
#define LONF2S lonf2_
#define LATG2S latg2_
#endif			/* MP */
#else			/* RSM */
#include <rscomf.rerun.h>
#include <rscommap.h>
#define LONF2 igrd12_
#define LATG2 jgrd12_
#define ILEN igrd1_
#define JLEN jgrd1_
#ifdef MP
#define MPGP2F rmpgp2f
#define MPABORT rmpabort
#include <comrmpi.h>
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else			/* MP */
#define MPABORT abort
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif			/* MP */
#endif			/* RSM */
c
c  downward sw fluxes from sw sib rad..
c 
      common/sibsw/ dfvbr(LONF2S,LATG2S),dfnbr(LONF2S,LATG2S),
     1              dfvdr(LONF2S,LATG2S),dfndr(LONF2S,LATG2S)
c
#ifdef CLR
      dimension cfsw(LONF2,LATG2,3)
      dimension cflw(LONF2,LATG2,3)
#endif
c
#ifdef DG3
#include <comgda.h>
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
      parameter(lonb=LONF2/2,latb=LATG2*2)
#ifdef DG3
      logical lbm(LONF2,LATG2)
      character g(200+LONF2*LATG2*(32+1)/8)
      dimension ids(255)
      dimension iens(5)
      dimension iclyr(3),ictl(3),icbl(3),itlcf(3)
      data iclyr/ihclyr,imclyr,ilclyr/
      data ictl /ihctl ,imctl ,ilctl /
      data icbl /ihcbl ,imcbl ,ilcbl /
      data itlcf/itoa,isfc,icolmn/
#ifdef MP
      real, allocatable :: fullat(:),fullon(:)
#endif
c
#ifndef RSM
      parameter(icen=7,icen2=0,igen=99)
#endif
c
#ifdef ASSIGN
      character*120 asgnstr
#endif
c
      dimension slmask(LONF2S,LATG2S)
c
      dimension work(LONF2,LATG2),slmsep(LONF2,LATG2),
     1          work2(LONF2,LATG2,levs_+2/levs_),
     1          workc(LONF2,20),iworkc(LONF2)
      dimension fluxw(LONF2,LATG2,26)
#ifdef MP
      dimension work2p(LONF2S,LATG2S,levs_+2/levs_)
#endif
c
      character*128 fno
#ifdef SCM
	  colrad(:)=90.
c#else
c      cl1=colrad(1)
#endif

c
#ifdef MP
      if(mype.eq.master) then
        print *,' start wridia '
#endif
#ifndef RSM
        call fnam('dia',3,thour,fno,ncho)
#else
        call fnam('r_dia',5,thour,fno,ncho)
#endif
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
      endif	! mype.eq.master
#endif
      call idsdef(1,ids)
      ienst=0
      iensi=0
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
      ifhr=nint(fhour)
      ithr=nint(thour)
c -- Never show 0-0hr averages in wgrib label, always round up to 1 hour
      if(ithr.eq.0) ithr=1
c
      dhour=dtpost/3600.
      if(dtpost.gt.0.) then
        rtime=1./dtpost
      else
        rtime=0.
      endif
c
#ifdef RSM
      colrad=0.
      proj=rproj
      delx=rdelx
      dely=rdely
#ifdef MP
      allocate (fullat(ILEN*JLEN),fullon(ILEN*JLEN))
      call MPGP2F(flat,igrd12p_,jgrd12p_,fullat,igrd12_,jgrd12_,1)
      call MPGP2F(flon,igrd12p_,jgrd12p_,fullon,igrd12_,jgrd12_,1)
#define FLAT fullat
#define FLON fullon
#else
#define FLAT flat
#define FLON flon
#endif
      call shalfo(FLAT,1)
      call shalfo(FLON,1)
      rlat1=FLAT(1)
      rlat2=FLAT(ILEN*JLEN)
      rlon1=FLON(1*1)
      rlon2=FLON(ILEN*JLEN)
#ifndef MP
      call shalfi(FLAT,1)
      call shalfi(FLON,1)
#else
      deallocate (fullat,fullon)
#endif
      nproj=rproj
      if( nproj.eq.0 ) then   ! mercater
        idrt=1
        ortru=rtruth
        proj=0.0
      elseif( abs(nproj).eq.1 ) then  ! polar projection
        idrt=5
        ortru=rorient
        proj=rproj
      else
        print *,' error in rwrtsfc for projection type....'
      endif
#else
#ifdef SCM
      idrt=0
#else
      idrt=4
#endif
      rlat1=0.
      rlon1=0.
      rlat2=0.
      rlon2=0.
      delx=0.
      dely=0.
      ortru=0.
      proj=0.
#endif
#ifdef MP
      call MPGP2F(slmask,LONF2S,LATG2S,slmsep,LONF2,LATG2,1)
      if(mype.eq.master) then
#else
      do j=1,LATG2
        do i=1,LONF2
          slmsep(i,j)=slmask(i,j)
        enddo
      enddo
#endif
#ifndef RSM
      call rowsep(slmsep)
#else
      call shalfo(slmsep,1)
#endif
#ifdef MP
      endif
#endif
c
#ifdef MP
      call MPGP2F(fluxr,LONF2S,LATG2S,fluxw,LONF2,LATG2,26)
      if(mype.eq.master) then
#else
      do k=1,26
       do j=1,LATG2
        do i=1,LONF2
          fluxw(i,j,k)=fluxr(i,j,k)
        enddo
       enddo
      enddo
#endif
#ifndef RSM
      do k=1,26
        call rowsep(fluxw(1,1,k))
      enddo
#else
      call shalfo(fluxw,26)
#endif
#ifdef CLR
      do j=1,LATG2
        do i=1,LONF2
          work(i,j)=fluxw(i,j,21)*rtime
        enddo
      enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'fluxw-21')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,iulwfc,itoa,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iulwfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
      do j=1,LATG2
        do i=1,LONF2
          work(i,j)=fluxw(i,j,22)*rtime
        enddo
      enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'fluxw-22')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,iuswfc,itoa,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
      do j=1,LATG2
        do i=1,LONF2
          work(i,j)=fluxw(i,j,25)*rtime
        enddo
      enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'fluxw-25')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idlwfc,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idlwfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
      do j=1,LATG2
        do i=1,LONF2
          work(i,j)=fluxw(i,j,23)*rtime
        enddo
      enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'fluxw-23')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idswfc,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
      do j=1,LATG2
        do i=1,LONF2
          work(i,j)=fluxw(i,j,24)*rtime
        enddo
      enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'fluxw-24')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,iuswfc,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(iuswfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
c   compute sw cloud forcing at toa (cld-clear),
c   flip sign so positive means cld is warming relative to clear
c
      do 600 j=1,LATG2
       do 600 i=1,LONF2
        cfsw(i,j,1) = - (fluxw(i,j,2)-fluxw(i,j,22))
  600 continue
c
c   compute cloud forcing at sfc (cld-clear)
c   again flip sign so positive means cld is warming relative to clear
c
      do 601 j=1,LATG2
       do 601 i=1,LONF2
       cfsw(i,j,2) =-(fluxw(i,j,3)-fluxw(i,j,4)
     1                      -(fluxw(i,j,24)-fluxw(i,j,23)))
  601 continue
c
c   flip sign so positive means cld is warming relative to clear
c
      do 602 j=1,LATG2
       do 602 i=1,LONF2
        cfsw(i,j,3) = - (cfsw(i,j,2) - cfsw(i,j,1))
  602 continue
c
      do 603 k=1,3
        do j=1,LATG2
          do i=1,LONF2
            work(i,j)=cfsw(i,j,k)*rtime
          enddo
        enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'cfsw')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,inswfc,itlcf(k),0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(inswfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
  603 continue
c
c   compute lw cloud forcing at toa (cld-clear),
c   flip sign so positive means cld is warming relative to clear
c
      do 1600 j=1,LATG2
       do 1600 i=1,LONF2
        cflw(i,j,1) = - (fluxw(i,j,1)-fluxw(i,j,21))
 1600 continue
c
c    compute cloud forcing at sfc (cld-clear)
c    again flip sign so positive means cld is warming relative to clear
c
      do 1601 j=1,LATG2
       do 1601 i=1,LONF2
       cflw(i,j,2) = - (fluxw(i,j,25)-fluxw(i,j,19))
 1601 continue
c
c   flip sign so positive means cld is warming relative to clear
c
      do 1602 j=1,LATG2
       do 1602 i=1,LONF2
        cflw(i,j,3) = - (cflw(i,j,2) - cflw(i,j,1))
 1602 continue
c
      do 1603 k=1,3
       do j=1,LATG2 
         do i=1,LONF2
           work(i,j)=cflw(i,j,k)
           work(i,j)=work(i,j)*rtime
         enddo
       enddo
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'cflw')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,inlwfc,itlcf(k),0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(inlwfc),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
 1603 continue
#endif			/* CLR */
c
#ifdef MP
      endif ! mype.eq.master
      call MPGP2F(dfvbr,LONF2S,LATG2S,work,LONF2,LATG2,1)
      if(mype.eq.master) then
#endif
      do j=1,LATG2
        do i=1,LONF2
#ifdef MP
          work(i,j)= work(i,j)*rtime
#else
          work(i,j)=dfvbr(i,j)*rtime
#endif
        enddo
      enddo
#ifndef RSM
      call rowsep(work)
#else
      call shalfo(work,1)
#endif
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'dfvbr')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idswvb,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswvb),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
#ifdef MP
      endif ! mype.eq.master
      call MPGP2F(dfvdr,LONF2S,LATG2S,work,LONF2,LATG2,1)
      if(mype.eq.master) then
#endif
      do j=1,LATG2
        do i=1,LONF2
#ifdef MP
          work(i,j)= work(i,j)*rtime
#else
          work(i,j)=dfvdr(i,j)*rtime
#endif
        enddo
      enddo
#ifndef RSM
      call rowsep(work)
#else
      call shalfo(work,1)
#endif
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'dfvdr')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idswvd,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswvd),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
#ifdef MP
      endif ! mype.eq.master
      call MPGP2F(dfnbr,LONF2S,LATG2S,work,LONF2,LATG2,1)
      if(mype.eq.master) then
#endif
      do j=1,LATG2
        do i=1,LONF2
#ifdef MP
          work(i,j)= work(i,j)*rtime
#else
          work(i,j)=dfnbr(i,j)*rtime
#endif
        enddo
      enddo
#ifndef RSM
      call rowsep(work)
#else
      call shalfo(work,1)
#endif
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'dfnbr')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idswnb,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswnb),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
c
#ifdef MP
      endif ! mype.eq.master
      call MPGP2F(dfndr,LONF2S,LATG2S,work,LONF2,LATG2,1)
      if(mype.eq.master) then
#endif
      do j=1,LATG2
        do i=1,LONF2
#ifdef MP
          work(i,j)= work(i,j)*rtime
#else
          work(i,j)=dfndr(i,j)*rtime
#endif
        enddo
      enddo
#ifndef RSM
      call rowsep(work)
#else
      call shalfo(work,1)
#endif
#ifdef DBG
      call maxmin(work,LONF2*LATG2,1,1,1,'dfndr')
#endif
      call gribit(work,lbm,idrt,lonb,latb,16,colrad,28,2,icen,igen,
     &            0,idswnd,isfc,0,0,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(idswnd),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)

#ifdef MP
      endif ! mype.eq.master
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
      call MPGP2F(work2p,LONF2S,LATG2S,work2,LONF2,LATG2,levs_)
      if(mype.eq.master) then
#endif
      ipu=ipugda(kd)
      ibm=ibmgda(kd)
      do 730 k=1,levs_
      isl=nint(sl(k)*1.e4)
#ifndef RSM
      call rowsep(work2(1,1,k))
#else
      call shalfo(work2(1,1,k),1)
#endif
#ifdef DBG
      print *,'kd=',kd,' k=',k
      call maxmin(work2(1,1,k),LONF2*LATG2,1,1,1,'diags')
#endif
      ibitmap=0
      if(ibm.ne.0) then
        ibitmap=1
        do j=1,LATG2
          do i=1,LONF2
            lbm(i,j)=work2(i,j,k).ne.0.
          enddo
        enddo
      endif
      call gribit(work2(1,1,k),lbm,idrt,lonb,latb,16,colrad,28,2,icen,
c    &            ibitmap,ipu,isglev,0,isl,iyr,imo,ida,ihr,
     &            igen,0,ipu,isglev,0,isl,iyr,imo,ida,ihr,
     &            iftime,ifhr,ithr,iavg,0,0,icen2,ids(ipu),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,
     &            g,lg,ierr)
      if(ierr.eq.0) call wryte(ndg,lg,g)
730   continue
c
#ifdef MP
      endif ! mype.eq.master
#endif
750   continue
c
#ifdef MP
      if(mype.eq.master) then
#endif
      close(ndg)
#ifdef MP
      endif
#endif
c
#endif
      return
      end
