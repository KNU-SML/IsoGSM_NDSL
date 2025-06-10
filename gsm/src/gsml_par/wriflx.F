#define IO2 io2
#define JO2 jo2
#define HVAP hvap
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
#define LNT2 lnt2
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
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#include <define.h>
#ifndef COUPLE_OUT
      subroutine wriflx(nn)
#else
      subroutine wriflx(nn,nnp)
#endif

#include <paramodel.h>
#ifndef RSM
#include <comio.h>
#endif

#ifdef RSM
#include <rscomf.rerun.h>
#include <rscommap.h>
#ifdef MP
#include <comrmpi.h>
#define MPABORT rmpabort
#else
#define MPABORT abort
#endif
#else
#include <comfphys.h>
#include <comfspec.vr.h>
#include <comfver.h>
#include <comfgrid.h>
#ifdef ISOTOPE
#include <comfiso.h>
#endif
#include <radiag.h>
#ifdef MP
#include <commpi.h>
#define MPABORT mpabort
#else
#define MPABORT abort
#endif
#endif
      parameter(iprs=1,itemp=11,iznlw=33,imerw=34,isphum=51,ipwat=54,
     $          ipcpr=59,isnowd=65,icldf=71,iccldf=72,
     $          islmsk=81,izorl=83,ialbdo=84,isoilm=144,icemsk=91,
     $          ilhflx=121,ishflx=122,izws=124,imws=125,ighflx=155,
     $          iuswfc=160,idswfc=161,iulwfc=162,idlwfc=163,
     $          inswfc=164,inlwfc=165,
     $          idswvb=166,idswvd=167,idswnb=168,idswnd=169,
     $          itmx=15,itmn=16,irnof=90,iep=145,
     $          isnwfll=64,isnwevp=230,isnwmlt=229,
     &          icldwk=146,izgw=147,imgw=148,ihpbl=221,
     $          idswf=204,idlwf=205,iuswf=211,iulwf=212,icpcpr=214,
     $          icanopy=223,iqull=202,iqvll=203,
     $          ivegtyp=225,ivegcov=87,isoltyp=154,ifwidx=209,
#ifdef ISOTOPE
     $          ipcpri1=67,ipcpri2=68,icpcpri1=69,icpcpri2=70,
     $          icondi1=112,icondi2=113,iccondi1=115,iccondi2=116,
     $          icond=111,iccond=114,
     $          ieqfo=157,ieqfd=158,
     $          ilhflxi1=119, ilhflxi2=120,irresv1=94,irresv2=95,
     $          ipwati1=150,ipwati2=151,isphumi1=153,isphumi2=152,
     $          irnofi1=92,irnofi2=93,ibgruni1=226,ibgruni2=227,
#endif
#ifdef CO2PROG
** co2 **
     $          ilhflxi1=119,
#endif
#ifdef RIVER
     $          igdriv=21,irflow=22,iimap=23,iroff=24,
#ifdef ISOTOPE
     $          igdrivi1=25,irflowi1=26,iroffi1=27,
     $          igdrivi2=28,irflowi2=29,iroffi2=30,
#endif
#endif
#ifndef NOAHYDRO
     $          ialhtfl=236,ievcnp=180,ibgrun=234)
#else
     $          ialhtfl=236,ievcnp=180,ibgrun=234,igeshem2=255)
#endif
      parameter(isfc=1,itoa=8,ielev=105,
     $          isglev=107,idbls=111,i2dbls=112,icolmn=200,
     $          ilcbl=212,ilctl=213,ilclyr=214,
     $          imcbl=222,imctl=223,imclyr=224,
     $          ihcbl=232,ihctl=233,ihclyr=234)
      parameter(inst=10,iavg=3,iacc=4)
      parameter(ifhour=1,ifday=2)
c
#ifdef MP
#ifdef RSM
#define MPGP2F rmpgp2f
#define ILEN2S igrd12p_
#define JLEN2S jgrd12p_
#define ILEN igrd1_
#define ILEN2 igrd12_
#define JLEN jgrd1_
#define JLEN2 jgrd12_
#else
#define MPGP2F mpgp2f
#define ILEN2S lonf2p_
#define JLEN2S latg2p_
#define ILEN lonf_
#define ILEN2 lonf2_
#define JLEN latg_
#define JLEN2 latg2_
#endif
#else
#ifdef RSM
#define ILEN2S igrd12_
#define JLEN2S jgrd12_
#define ILEN igrd1_
#define ILEN2 igrd12_
#define JLEN jgrd1_
#define JLEN2 jgrd12_
#else
#define ILEN2S lonf2_
#define JLEN2S latg2_
#define ILEN lonf_
#define ILEN2 lonf2_
#define JLEN latg_
#define JLEN2 latg2_
#endif
#endif
c
c  full arrays
c
#ifdef NOAHYDRO
      parameter(hvap=hvap_)
#endif
      parameter(len=ILEN*JLEN)
      logical   lbm(len)
      character g(200+len*(32+1)/8)
c
      parameter(nfld=16)
      integer   ipur(nfld),itlr(nfld)
      data ipur/iulwf , iuswf , iuswf , idswf ,  icldf,   iprs,
     $          iprs, itemp ,  icldf,   iprs,   iprs, itemp ,
     $          icldf,   iprs,   iprs, itemp /
      data itlr/itoa  , itoa  , isfc  , isfc  , ihclyr, ihctl ,
     $          ihcbl , ihctl , imclyr, imctl , imcbl , imctl ,
     $          ilclyr, ilctl , ilcbl , ilctl /
      real rtimer(nfld)
c
      parameter(lens=ILEN2S*JLEN2S)
c
c  full working arrays
c
      dimension work(len,4)
      dimension slmsep(len)
#ifdef MP
      dimension fluxf(lens,4)
      dimension workp(lens)
#else
      dimension fluxf(len,4)
#endif
c
      dimension lsoimlev(lsoil_+1)
#ifdef OSULSM1
      data lsoimlev/0,10,200/
#endif
#ifdef OSULSM2
      data lsoimlev/0,10,200/
#endif
#ifdef NOALSM1
      data lsoimlev/0,10,40,100,200/
#endif
#ifdef VIC
      data lsoimlev/0,1,2,3/
#endif
#ifdef VIC
      dimension lsoitlev(nsoil_+1)
      data lsoitlev/0,1,2,3,4,5/
#endif
c
#ifdef MP
      real, allocatable :: fullat(:),fullon(:)
#endif
c
      dimension ids(255)
      dimension iens(5)
#ifndef RSM
#undef MRG_POST
#endif
c
      character*128 fno
#ifdef ASSIGN
      character*128 asgnstr
#endif
c
#ifdef RIVER
#include <comfriv.h>
c
      parameter(lenr=io2_*(jo2_))
      logical   lbmr(lenr)
      character gr(200+lenr*(16+1)/8)
      real gdriv0(io2_,jo2_)
      real rflow0(io2_,jo2_)
#endif
c
#ifdef RAD_SMOOTH_CLOUD
      print *,'this version does not work for rad_smooth_cloud'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
#endif
c
      if(nn.le.0) return
c
#ifdef MP
      if(mype.eq.master) then
#endif
#ifdef RSM
#ifndef MRG_POST
        call fnam('r_flx',5,thour,fno,ncho)
#endif
#else
        call fnam('flx',3,thour,fno,ncho)
#endif
#ifndef MRG_POST
#ifdef ASSIGN
        write(asgnstr,'(23hassign -s unblocked  u:,I2,)') nn
        call assign('assign -R')
        call assign(asgnstr)
#endif
        open(unit=nn,file=fno(1:ncho),form='unformatted',err=900)
        go to 901
  900   continue
        write(6,*) 'wriflx: error in opening file ',fno(1:ncho)
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
  901   continue
#ifndef NOPRINT
        write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
        rewind nn
#endif
#ifdef MP
      endif
#endif
c
      call idsdef(1,ids)
      ilpds=28
      if(icen2.eq.2) ilpds=45
c
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
      ifhr=nint(fhour)
      ithr=nint(thour)
c -- Never show 0-0hr averages, always round up to 1 hour
      if(ithr.eq.0) ither=1
c
      dhour=dtpost/3600.
      if(dtpost.gt.0) then
        rtime=1./dtpost
      else
        rtime=0.
      endif
c
      secswr=max(dhour,dtswav) * 3600.
      seclwr=max(dhour,dtlwav) * 3600.
      if(secswr.gt.0.) then
        rtimsw=1./secswr
      else
        rtimsw=1.
      endif
      if(seclwr.gt.0.) then
        rtimlw=1./seclwr
      else
        rtimlw=1.
      endif
      do n=1,nfld
      rtimer(n)=rtimsw
      enddo
      rtimer(1)=rtimlw
c
#ifdef RSM
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
c
#ifdef MP
      call MPGP2F(slmsk,ILEN2S,JLEN2S,slmsep,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      slmsep(n)=slmsk(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(slmsep,1)
#else
      call rowsep(slmsep)
#ifdef SCM
      cl1=90.
#else
      cl1=colrad(1)
#endif
#endif
#ifdef MP
      endif
#endif
c
#ifndef RSM
#define MRG_POST
#endif
#ifndef MRG_POST
#define GRIBIT  write(nn)
#define GRIBIT2  write(nn)
#else
#define GRIBIT  call gribit(
#define GRIBIT2  call gribit2(
#endif
c
c 1  surface u-stress
c
#ifdef MP
      call MPGP2F(dusfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dusfc(n,1)*rtime
#endif
      enddo
      maxbit=16
#ifndef NOAHYDRO
c      maxbit2=16
      maxbit2=32
#else
      maxbit2=32
#endif
      colat=0.
      iptv=2
      ibm0=0
      ibm1=1
      il1k=0
      il2k=0
      ip1=0
      ip2=0
      ina=0
      inm=0
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'uflx')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,izws,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(izws),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'uflx grib created    '
#endif
      else
        print *,'uflx grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 2  surface v-stress
c
#ifdef MP
      call MPGP2F(dvsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dvsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'vflx')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,imws,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(imws),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj 
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'vflx grib created    '
#endif
      else
        print *,'vflx grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 3  sensible heat flux
c
#ifdef MP
      call MPGP2F(dtsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dtsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'tflx')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ishflx,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ishflx),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'tflx grib created    '
#endif
      else
        print *,'tflx grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 4  latent heat flux
c
#ifdef MP
      call MPGP2F(dqsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dqsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflx')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ilhflx,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ilhflx),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'qflx grib created    '
#endif
      else
        print *,'qflx grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 5  surface temperature
c
#ifdef MP
      call MPGP2F(tsea,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=tsea(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'tsfc')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itemp,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itemp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'tsfc grib created    '
#endif
      else
        print *,'tsfc grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 6  soil moisture (multi-level)
c
      do k=1,lsoil_
#ifdef MP
        call MPGP2F(smc(1,1,k),ILEN2S,JLEN2S,work,
     1              ILEN2,JLEN2,1)
        if(mype.eq.master) then
#else
        do n=1,len
          work(n,1)=smc(n,1,k)
        enddo
#endif
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
        il1k=lsoimlev(k)
        il2k=lsoimlev(k+1)
#ifdef DBG
        print *,'soil level=',k
        call maxmin(work,ILEN*JLEN,1,1,1,'soilw')
#endif
        GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isoilm,i2dbls,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(isoilm),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'soilm level',k,' grib created'
#endif
        else
          print *,'soilm level',k,' grib make failed'
        endif
#endif
#ifdef MP
        endif
#endif
      enddo
      il1k=0
      il2k=0
c
c 8  soil temperature (multi-level)
c
#ifndef VIC
      do k=1,lsoil_
#else
      do k=1,nsoil_
#endif
#ifdef MP
        call MPGP2F(stc(1,1,k),ILEN2S,JLEN2S,work,
     1              ILEN2,JLEN2,1)
        if(mype.eq.master) then
#else
        do n=1,len
          work(n,1)=stc(n,1,k)
        enddo
#endif
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef VIC
        il1k=lsoitlev(k)
        il2k=lsoitlev(k+1)
#else
        il1k=lsoimlev(k)
        il2k=lsoimlev(k+1)
#endif
#ifdef DBG
        print *,'soil level=',k
        call maxmin(work,ILEN*JLEN,1,1,1,'soilt')
#endif
        GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itemp,i2dbls,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itemp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'soil temperature level 1 grib created    '
#endif
        else
          print *,'soil temperature level 1 grib make failed'
        endif
#endif
#ifdef MP
        endif
#endif
      enddo
      il1k=0
      il2k=0
c
c 10  soil water equivalent snow depth
c
#ifdef MP
      call MPGP2F(sheleg,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=sheleg(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'snow')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isnowd,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(isnowd),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'snow depth grib created    '
#endif
      else
        print *,'snow depth grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 11 downward longwave radiation flux at the surface
c
#ifdef MP
      call MPGP2F(dlwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dlwsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'dlwsfc')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,idlwf,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(idlwf),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'dlwflx sfc grib created    '
#endif
      else
        print *,'dlwflx sfc grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 12 upword longwave radiation flux at the surface
c
#ifdef MP
      call MPGP2F(ulwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=ulwsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ulwsfc')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iulwf,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iulwf),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'ulwflx sfc grib created    '
#endif
      else
        print *,'ulwflx sfc grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif

c
c   13 upward long wave radiation flux at the top of the atmosphere
c   14 upward short wave radiation flux at the top of the atmosphere
c   15 upward short wave radiation flux at the surface
c   16 downward short wave radiation flux at the surface
c
#ifdef MP
#ifdef RSMVECTORIZE
      call fillarray(fluxr,26)
#endif
#endif
      do 113 k=1,4
        ii=0
#ifdef MP
        call MPGP2F(fluxr(1,1,k),ILEN2S,JLEN2S,work,
     1              ILEN2,JLEN2,1)
        if(mype.eq.master) then
#else
        do j=1,JLEN2
          do i=1,ILEN2
            ii=ii+1
            work(ii,1)=fluxr(i,j,k)
          enddo
        enddo
#endif
        do n=1,len
          work(n,1)=work(n,1)*rtimer(k)
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'radflxes')
#endif
        GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipur(k),itlr(k),il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(k)),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &            ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'radiation fluxes grib created    '
#endif
        else
          print *,'radiation fluxes grib make failed'
        endif
#endif
#ifdef COUPLE_OUT
        if ( k.eq.3.or.k.eq.4) then
          write(nnp) (work(n,1),n=1,len)
        endif
#endif
#ifdef MP
         endif
#endif
  113   continue
c
      iw1=1
      iw2=iw1+4*ILEN2
      iw3=iw2+4*ILEN2
      iw4=iw3+4*ILEN2
      iw5=iw4+4*ILEN2
      do 813 k=5,7
        ij=0
        do j=1,JLEN2S
         do i=1,ILEN2S
          ij=ij+1
          if(fluxr(i,j,k).gt.0.) then
            fluxf(ij,2) = fluxr(i,j,k+3) / fluxr(i,j,k)
            fluxf(ij,3) = fluxr(i,j,k+6) / fluxr(i,j,k)
            fluxf(ij,4) = fluxr(i,j,k+9) / fluxr(i,j,k)
            fluxf(ij,1) = fluxr(i,j,k) * rtimsw
          else
c  zero cld top temp if no clds--safety, cause use zero in ggintt
            fluxr(i,j,k+9) = 0.
            fluxf(ij,2) = fluxr(i,j,k+3)
            fluxf(ij,3) = fluxr(i,j,k+6)
            fluxf(ij,4) = fluxr(i,j,k+9)
            fluxf(ij,1) = fluxr(i,j,k) * rtimsw
          endif
         enddo
        enddo
#ifdef MP
        call MPGP2F(fluxf,ILEN2S,JLEN2S,work,ILEN2,JLEN2,4)
        if(mype.eq.master) then
#endif
c
#ifdef MP
#define WORK work
#else
#define WORK fluxf
#endif

c
c   17 high level cloud amount
c   21 middle level cloud amount
c   25 low level cloud amount
c
        k4=4+(k-5)*4
#ifdef RSM
        call shalfo(WORK(1,1),1)
#else
        call rowsep(WORK(1,1))
#endif
        do n=1,len
          WORK(n,1)=WORK(n,1)*1.e2
        enddo
        l=k4+1
#ifdef DBG
        call maxmin(WORK,ILEN*JLEN,1,1,1,'cloud amount')
#endif
        GRIBIT WORK(1,1),lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG

          print *,'cloud amount grib created    '
#endif
        else
          print *,'cloud amount grib make failed'
        endif
#endif
c
c   18 high level cloud top pressure
c   22 middle level cloud top pressure
c   26 low level cloud top pressure
c
#ifdef RSM
        call shalfo(WORK(1,2),1)
#else
        call rowsep(WORK(1,2))
#endif
        do n=1,len
          WORK(n,2)=WORK(n,2)*1.e3
        enddo
        l=k4+2
#ifdef DBG
        call maxmin(WORK(1,2),ILEN*JLEN,1,1,1,'cloud top pres')
#endif
        GRIBIT WORK(1,2),lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'cloud top pressure grib created    '
#endif
        else
          print *,'cloud top pressure grib make failed'
        endif
#endif
c
c   19 high level cloud bottom pressure
c   23 middle level cloud bottom pressure
c   27 low level cloud bottom pressure
c
#ifdef RSM
        call shalfo(WORK(1,3),1)
#else
        call rowsep(WORK(1,3))
#endif
        do n=1,len
          WORK(n,3)=WORK(n,3)*1.e3
        enddo
        l=k4+3
#ifdef DBG
        call maxmin(WORK(1,3),ILEN*JLEN,1,1,1,'cloud bot pres')
#endif
        GRIBIT WORK(1,3),lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'cloud bottom pressure grib created    '
#endif
        else
          print *,'cloud bottom pressure grib make failed'
        endif
#endif
c
c   20 high level cloud top tempeature
c   24 middle level cloud top tempeature
c   28 low level cloud top tempeature
c
#ifdef RSM
        call shalfo(WORK(1,4),1)
#else
        call rowsep(WORK(1,4))
#endif
        l=k4+4
#ifdef DBG
        call maxmin(WORK(1,4),ILEN*JLEN,1,1,1,'cloud top tmp')
#endif
        GRIBIT WORK(1,4),lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipur(l),itlr(l),il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipur(l)),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
        if(ierr.eq.0) then
          call wryte(nn,lg,g)
#ifdef DBG
          print *,'cloud top temp grib created    '
#endif
        else
          print *,'cloud top temp grib make failed'
        endif
#endif
c
#ifdef MP
        endif
#endif
  813 continue
c
c 29  precipitation rate
c
#ifdef MP
      call MPGP2F(geshem,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshem(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'prate')
#endif
#ifndef ES
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
#else
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
#endif
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipcpr,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipcpr),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'precip grib created    '
#endif
      else
        print *,'precip grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 30  convective precpipitation rate
c
#ifdef MP
      call MPGP2F(bengsh,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengsh(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cprate')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icpcpr,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icpcpr),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective precip grib created    '
#endif
      else
        print *,'convective precip grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 31  ground flux
c
#ifdef MP
      call MPGP2F(gflux,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=gflux(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'gflux')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ighflx,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ighflx),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'ground heat flux grib created    '
#endif
      else
        print *,'ground heat flux grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 32  land sea mask
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do n=1,len
        work(n,1)=mod(slmsep(n),2.)
      enddo
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'lsmask')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,islmsk,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(islmsk),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'land sea mask grib created    '
#endif
      else
        print *,'land sea mask grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 33  ice mask
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do n=1,len
        work(n,1)=max(slmsep(n)-1.,0.)
      enddo
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ice mask')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icemsk,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(icemsk),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'ice mask grib created    '
#endif
      else
        print *,'ice mask grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 34  10m wind u
c
#ifdef MP
      call MPGP2F(u10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=u10m(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=10
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'u10m')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iznlw,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(iznlw),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'10m u grib created    '
#endif
      else
        print *,'10m u grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 35  10m wind v
c
#ifdef MP
      call MPGP2F(v10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=v10m(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=10
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'v10m')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,imerw,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(imerw),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'10m v grib created    '
#endif
      else
        print *,'10m v grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 36  temperature 2m
c
#ifdef MP
      call MPGP2F(t2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=t2m(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'t2m')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itemp,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itemp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'2m temp grib created    '
#endif
      else
        print *,'2m temp grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 37  specific humidity 2m
c
#ifdef MP
      call MPGP2F(q2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=q2m(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'q2m')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isphum,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(isphum),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'2m specific humidity grib created    '
#endif
      else
        print *,'2m specific humidity grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 39  max temperature
c
#ifdef MP
      call MPGP2F(tmpmax,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=tmpmax(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'max temp')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itmx,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itmx),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'max temp grib created    '
#endif
      else
        print *,'max temp grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c  reset tmpmax
c
      do n=1,lens
      tmpmax(n,1) = 0.
      enddo
c
c 40  min temperature
c
#ifdef MP
      call MPGP2F(tmpmin,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=tmpmin(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'min temp')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itmn,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itmn),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'min temp grib created    '
#endif
      else
        print *,'min temp grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c  reset tmpmin
c
      do n=1,lens
      tmpmin(n,1) = 1.e10
      enddo
c
c 41  runoff
c
#ifdef MP
      call MPGP2F(runoff,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=runoff(n,1) * 1.e3 * rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'runoff')
#endif
#ifndef ES
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
#else
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
#endif
     &            ilpds,iptv,icen,igen,
     &            ibm0,irnof,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irnof),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'runoff grib created    '
#endif
      else
        print *,'runoff grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 42  potential evaporation
c
#ifdef MP
      call MPGP2F(ep,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=ep(n,1) * rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'potevap')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iep,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iep),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'potential evap grib created    '
#endif
      else
        print *,'potential evap grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
#ifdef SAS
c
c 43  cloud work function
c
c#ifdef MP
c      call MPGP2F(cldwrk,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
c      if(mype.eq.master) then
c#endif
c      do n=1,len
c#ifdef MP
c      work(n,1)=work(n,1)*rtime
c#else
c      work(n,1)=cldwrk(n,1) * rtime
c#endif
c      enddo
c#ifdef RSM
c      call shalfo(work,1)
c#else
c      call rowsep(work)
c#endif
c#ifdef DBG
c      call maxmin(work,ILEN*JLEN,1,1,1,'cloud wk func.')
c#endif
c      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
c     &            ilpds,iptv,icen,igen,
c     &            ibm0,icldwk,isfc,il1k,il2k,iyr,imo,ida,ihr,
c     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icldwk),iens,
c     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
c#ifdef MRG_POST
c     &           ,g,lg,ierr)
c      if(ierr.eq.0) then
c        call wryte(nn,lg,g)
c#ifdef DBG
c        print *,'cloud work func grib created    '
c#endif
c      else
c        print *,'cloud work func grib make failed'
c      endif
c#endif
c#ifdef MP
c      endif
c#endif
#endif
c
c 44  gravity wave drag (u)
c
#ifdef MP
      call MPGP2F(dugwd,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dugwd(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ugwd')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,izgw,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(izgw),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'gravity wave drag u grib created    '
#endif
      else
        print *,'gravity wave drag u grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 45  gravity wave drag (v)
c
#ifdef MP
      call MPGP2F(dvgwd,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dvgwd(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'vgwd')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,imgw,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(imgw),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'gravity wave drag v grib created    '
#endif
      else
        print *,'gravity wave drag v grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 46  pbl height
c
#ifdef MP
      call MPGP2F(hpbl,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=hpbl(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pblh')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ihpbl,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ihpbl),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'pbl height grib created    '
#endif
      else
        print *,'pbl height grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 47  precipitable water
c
#ifdef MP
      call MPGP2F(pwat,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=pwat(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pwater')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipwat,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ipwat),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'precipitable water grib created    '
#endif
      else
        print *,'precipitable water grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 48  snow fall amount
c
#ifdef MP
      call MPGP2F(snowfall,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*1.e3*rtime
#else
       work(n,1)=snowfall(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'snowfall')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isnwfll,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwfll),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'snow fall grib created    '
#endif
      else
        print *,'snow fall grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 49  snow evaporation
c
#ifdef MP
      call MPGP2F(snowevap,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=snowevap(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'snowevap')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isnwevp,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwevp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'snow evap grib created    '
#endif
      else
        print *,'snow evap grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 50  snow melt
c
#ifdef MP
      call MPGP2F(snowmelt,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=snowmelt(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'snow melt')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isnwmlt,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isnwmlt),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'snow melt grib created    '
#endif
      else
        print *,'snow melt grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 51  moisture flux (u)
c
#ifdef MP
      call MPGP2F(qull,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=qull(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qfulx-u')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iqull,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iqull),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'moisture flux u grib created    '
#endif
      else
        print *,'moisture flux u grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 52  moisture flux (v)
c
#ifdef MP
      call MPGP2F(qvll,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=qvll(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflux-v')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iqvll,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iqvll),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'moisture flux v grib created    '
#endif
      else
        print *,'moisture flux v grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 53  canopy water content
c
#ifdef MP
      call MPGP2F(canopy,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)
#else
       work(n,1)=canopy(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'canopy water')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icanopy,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(icanopy),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'canopy water content grib created    '
#endif
      else
        print *,'canopy water content grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 54 downward short wave flux at the top of the atmosphere
c
#ifdef MP
      call MPGP2F(fluxr(1,1,18),ILEN2S,JLEN2S,work,
     1            ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      ii=0
      do j=1,JLEN2
        do i=1,ILEN2
          ii=ii+1
          work(ii,1)=fluxr(i,j,18)
        enddo
      enddo
#endif
      do n=1,len
        work(n,1)=work(n,1)*rtime
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'dswtoa')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,idswf,itoa,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(idswf),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'dswflx top grib created    '
#endif
      else
        print *,'dswflx top grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 55  total cloud cover
c
#ifdef MP
      call MPGP2F(fluxr(1,1,26),ILEN2S,JLEN2S,work,
     1            ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      ii=0
      do j=1,JLEN2
        do i=1,ILEN2
          ii=ii+1
          work(ii,1)=fluxr(i,j,26)
        enddo
      enddo
#endif
#ifdef RAD_SMOOTH_CLOUD
      do ij=1,ILEN2*JLEN2
        work(ij,1) = work(ij,1) * rtimsw
      enddo
      call cvinamt(work,
     1              ILEN2,JLEN2,JLEN2,
     2              fluxf(1,1),
     3              ILEN2,JLEN2,JLEN2,
     5              workc(iw1),workc(iw2),
     6              workc(iw5),iworkc,1,1,1)
c
c   reset data
c
      do ij=1,ILEN2*JLEN2
        work(ij,1) = work(ij,1) / rtimsw
      enddo
#else
      do ij=1,ILEN2*JLEN2
        work(ij,1) = work(ij,1) * rtimsw
      enddo
#endif
c
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work,1)
#endif
      do n=1,len
        work(n,1)=work(n,1)*1.e2
      enddo
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'total cloud')
#endif
      GRIBIT work(1,1),lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icldf,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icldf),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'total cloud cover grib created    '
#endif
      else
        print *,'total cloud cover grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 56  albedo
c
#ifdef MP
      call MPGP2F(fluxr(1,1,17),ILEN2S,JLEN2S,work,
     1            ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      ii=0
      do j=1,JLEN2S
        do i=1,ILEN2S
          ii=ii+1
          work(ii,1)=fluxr(i,j,17)
        enddo
      enddo
#endif
      do n=1,len
        work(n,1)=work(n,1)*rtimsw * 100.
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'albedo')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ialbdo,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ialbdo),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'albedo grib created    '
#endif
      else
        print *,'albedo grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 57  surface roughness
c
#ifdef MP
      call MPGP2F(zorl,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)
#else
       work(n,1)=zorl(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'z0')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,izorl,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(izorl),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'albedo grib created    '
#endif
      else
        print *,'surface roughness grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 58  vegetation type
c
#ifndef OSULSM1
#ifdef MP
      call MPGP2F(vtype,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)
#else
       work(n,1)=vtype(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'vegtype')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ivegtyp,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ivegtyp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'vegitation type grib created    '
#endif
      else
        print *,'vegitation type grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
#endif
c
c 59  vegetation cover
c
#ifdef OSULSM1
#define VARIAB plantr
#else
#define VARIAB vfrac
#endif

#ifdef MP
      call MPGP2F(VARIAB,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)
#else
       work(n,1)=VARIAB(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'vegcover')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ivegcov,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ivegcov),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'vegetation cover grib created    '
#endif
      else
        print *,'vegetation cover grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
#ifndef VIC
#ifndef OSULSM1
c
c 60  soil type
c
#ifdef MP
      call MPGP2F(stype,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)
#else
       work(n,1)=stype(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'soiltype')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isoltyp,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(isoltyp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'soil type grib created    '
#endif
      else
        print *,'soil type grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
#endif   /*ifndef OSULSM */
#endif   /* ifdef VIC */
c
c 61  Another latent heat flux (from progtm99)
c
#ifdef MP
      call MPGP2F(alhtfl,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=alhtfl(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflx progtm')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ialhtfl,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ialhtfl),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'another latent heat flux grib created    '
#endif
      else
        print *,'another latent heat flux grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 62  Canopy evaporation
c
#ifdef HYDRO
#define EVCNP evcnp
#else
#ifndef NOAHYDRO
#define EVCNP ep
#else
#define EVCNP evcnp
#endif
#endif
#ifdef MP
      call MPGP2F(EVCNP,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
#ifndef NOAHYDRO
       work(n,1)=work(n,1)*rtime
#else
       work(n,1)=work(n,1)*rtime*1.e3*hvap
#endif
#else
#ifndef NOAHYDRO
       work(n,1)=EVCNP(n,1)*rtime
#else
       work(n,1)=EVCNP(n,1)*rtime*1.e3*hvap
#endif
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'canopy evap')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ievcnp,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ievcnp),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'canopy evap grib created    '
#endif
      else
        print *,'canopy evap grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 62  Baseflow Groundwater Runoff
c
#ifdef MP
      call MPGP2F(bgrun,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*1.e03*rtime
#else
       work(n,1)=bgrun(n,1)*1.e03*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'base flow')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ibgrun,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ibgrun),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'bgrun grib created    '
#endif
      else
        print *,'bgrun grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 63  Surface Pressure
c
#ifdef MP
      call MPGP2F(psurf,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3
#else
      work(n,1)=psurf(n,1)*1.e3
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'prate')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iprs,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(iprs),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'surface pressure grib created    '
#endif
      else
        print *,'surface pressure grib make failed'
      endif
#endif
#ifdef COUPLE_OUT
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
#ifdef NOAHYDRO
c
c 64  Second Precip
c
#ifdef MP
      call MPGP2F(geshem2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*1.e03*rtime
#else
       work(n,1)=geshem2(n,1)*1.e03*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'base flow')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,igeshem2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(igeshem2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'geshem2 grib created    '
#endif
      else
        print *,'geshem2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
#endif
c
c 65  fire weather index
c
#ifndef SCM
#ifdef MP
#define WORKP workp
#else
#define WORKP work
#endif
      call fwindx(ILEN2S,JLEN2S,u10m,v10m,t2m,q2m,psurf,WORKP)
#ifdef MP
      call MPGP2F(WORKP,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'fwindx')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ifwidx,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ifwidx),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'fwindx grib created    '
#endif
      else
        print *,'fwindx grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c 
#ifndef NOPRINT
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'grib flux file created.'
#ifdef MP
      endif
#endif
c /*undef SCM */
#endif 
#ifdef CO2PROG
** co2 **
c
c 66 co2 flux 
c
#ifdef MP
      call MPGP2F(dqsfci1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dqsfci1(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflx1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ilhflxi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ilhflxi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'qflx1 grib created    '
#endif
      else
        print *,'qflx1 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c /* end of CO2PROG */
#endif
c
#ifdef ISOTOPE
c
c 66 isotope precipitation rate 1
c
#ifdef MP
      call MPGP2F(geshemi1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshemi1(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pratei1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipcpri1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipcpri1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'precipi1 grib created    '
#endif
      else
        print *,'precipi1 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 67  isotope precipitation rate 2
c
#ifdef MP
      call MPGP2F(geshemi2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshemi2(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pratei2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipcpri2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipcpri2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'precipi2 grib created    '
#endif
      else
        print *,'precipi2 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 68  convective precpipitation rate 1
c
#ifdef MP
      call MPGP2F(bengshi1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengshi1(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cpratei1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &          ilpds,iptv,icen,igen,
     &          ibm0,icpcpri1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icpcpri1),iens,
     &          rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective precipi1 grib created    '
#endif
      else
        print *,'convective precipi1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 69  convective precpipitation rate 2
c
#ifdef MP
      call MPGP2F(bengshi2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengshi2(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cpratei2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &          ilpds,iptv,icen,igen,
     &          ibm0,icpcpri2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icpcpri2),iens,
     &          rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective precipi2 grib created    '
#endif
      else
        print *,'convective precipi2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 70  latent heat flux isotope 1
c
#ifdef MP
      call MPGP2F(dqsfci1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dqsfci1(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflx1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ilhflxi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ilhflxi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'qflx1 grib created    '
#endif
      else
        print *,'qflx1 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 71  latent heat flux isotope 2
c
#ifdef MP
      call MPGP2F(dqsfci2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*rtime
#else
      work(n,1)=dqsfci2(n,1)*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'qflx2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ilhflxi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ilhflxi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'qflx2 grib created    '
#endif
      else
        print *,'qflx2 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 72  surface reservoir isotope ratio 1
c
#ifdef MP
      call MPGP2F(rresv1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=rresv1(n,1,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rresv1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irresv1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(irresv1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'rresv1 grib created    '
#endif
      else
        print *,'rresv1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 73  surface reservoir isotope ratio 2
c
#ifdef MP
      call MPGP2F(rresv2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=rresv2(n,1,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rresv2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irresv2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(irresv2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'rresv2 grib created    '
#endif
      else
        print *,'rresv2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 74  precipitable water isotope 1
c
#ifdef MP
      call MPGP2F(pwati1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=pwati1(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pwati1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipwati1,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ipwati1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'pwati1 grib created    '
#endif
      else
        print *,'pwati1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 75  precipitable water isotope 2
c
#ifdef MP
      call MPGP2F(pwati2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      work(n,1)=pwati2(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'pwati2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipwati2,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ipwati2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'pwati2 grib created    '
#endif
      else
        print *,'pwati2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 76  isotope specific humidity 2m
c
#ifdef MP
      call MPGP2F(q2mi1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=q2mi1(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'q2mi1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isphumi1,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(isphumi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'2m specific humidity grib created    '
#endif
      else
        print *,'2m specific humidity grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 77  isotope specific humidity 2m
c
#ifdef MP
      call MPGP2F(q2mi2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)
#else
      work(n,1)=q2mi2(n,1)
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
      il2k=2
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'q2mi2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,isphumi2,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(isphumi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'2m specific humidity grib created    '
#endif
      else
        print *,'2m specific humidity grib make failed'
      endif
#endif
      il2k=0
#ifdef MP
      endif
#endif
c
c 78  convective total condensation rate
c
#ifdef MP
      call MPGP2F(bengshc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengshc(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ccrate')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iccond,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iccond),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective cond grib created    '
#endif
      else
        print *,'convective cond grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 79  convective condensation rate 1
c
#ifdef MP
      call MPGP2F(bengshci1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengshci1(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ccratei1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &          ilpds,iptv,icen,igen,
     &          ibm0,iccondi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iccondi1),iens,
     &          rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective condi1 grib created    '
#endif
      else
        print *,'convective condi1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 80  convective condensation rate 2
c
#ifdef MP
      call MPGP2F(bengshci2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=bengshci2(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ccratei2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
     &          ilpds,iptv,icen,igen,
     &          ibm0,iccondi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &          ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iccondi2),iens,
     &          rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'convective condi2 grib created    '
#endif
      else
        print *,'convective condi2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 81  condensation rate
c
#ifdef MP
      call MPGP2F(geshemc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshemc(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'crate')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icond,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icond),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'cond grib created    '
#endif
      else
        print *,'cond grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 82  isotope condensation rate 1
c
#ifdef MP
      call MPGP2F(geshemci1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshemci1(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cratei1')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icondi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icondi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'condi1 grib created    '
#endif
      else
        print *,'condi1 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 83  isotope condensation rate 2
c
#ifdef MP
      call MPGP2F(geshemci2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=geshemci2(n,1)*1.e3*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cratei2')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icondi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(icondi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'condi2 grib created    '
#endif
      else
        print *,'condi2 grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 84  equilibrium fraction for 18o
c
#ifdef MP
      call MPGP2F(meqfo,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e2
#else
      work(n,1)=meqfo(n,1)*1.e2
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'meqfo')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ieqfo,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ieqfo),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'eqfo grib created    '
#endif
      else
        print *,'eqfo grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 85  equilibrium fraction for d
c
#ifdef MP
      call MPGP2F(meqfd,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e2
#else
      work(n,1)=meqfd(n,1)*1.e2
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'meqfd')
#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ieqfd,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ieqfd),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'eqfd grib created    '
#endif
      else
        print *,'eqfd grib make failed'
      endif
#endif
#ifdef COUPLE
      write(nnp) (work(n,1),n=1,len)
#endif
#ifdef MP
      endif
#endif
c
c 86  surface runoff isotope 1
c
#ifdef MP
      call MPGP2F(runoff1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=runoff1(n,1) * 1.e3 * rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'runoff1')
#endif
#ifndef ES
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
#else
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
#endif
     &            ilpds,iptv,icen,igen,
     &            ibm0,irnofi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irnofi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'runoff1 grib created    '
#endif
      else
        print *,'runoff1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 87  runoff isotope 2
c
#ifdef MP
      call MPGP2F(runoff2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
      work(n,1)=work(n,1)*1.e3*rtime
#else
      work(n,1)=runoff2(n,1) * 1.e3 * rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'runoff2')
#endif
#ifndef ES
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit2,colat,
#else
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
#endif
     &            ilpds,iptv,icen,igen,
     &            ibm0,irnofi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irnofi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'runoff grib created    '
#endif
      else
        print *,'runoff grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 88  Baseflow Groundwater Runoff isotope 1
c
#ifdef MP
      call MPGP2F(bgrun1,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*1.e03*rtime
#else
       work(n,1)=bgrun1(n,1)*1.e03*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
c#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'base flow1')
c#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ibgruni1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &           ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ibgruni1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'bgrun1 grib created    '
#endif
      else
        print *,'bgrun1 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
c
c 89  Baseflow Groundwater Runoff isotope 2
c
#ifdef MP
      call MPGP2F(bgrun2,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
      do n=1,len
#ifdef MP
       work(n,1)=work(n,1)*1.e03*rtime
#else
       work(n,1)=bgrun2(n,1)*1.e03*rtime
#endif
      enddo
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
c#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'base flow2')
c#endif
      GRIBIT work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ibgruni2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &           ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ibgruni2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
#ifdef MRG_POST
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'bgrun2 grib created    '
#endif
      else
        print *,'bgrun2 grib make failed'
      endif
#endif
#ifdef MP
      endif
#endif
*** Finish isotopes ***
#endif

#endif
c
#ifdef MP
      if(mype.eq.master) then
#endif
      close(nn)
#ifdef MP
      endif
#endif
c
#ifndef RSM
#undef MRG_POST
#endif
#ifdef RIVER
c
c kei's river
c
#ifndef RSM
#ifdef MP
      if(mype.eq.master) then
#endif
        call fnam('riv',3,thour,fno,ncho)
        open(unit=nn,file=fno(1:ncho),form='unformatted',err=902)
        go to 903
 902    continue
        write(6,*) ' error in opening file ',fno(1:ncho)
#ifdef MP
        call mpabort
#else
        call abort
#endif
 903    continue
#ifndef NOPRINT
        write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
        rewind nn
#ifdef MP
      endif
#endif
c
c 1 riv direction
c
#ifdef MP
      if(mype.eq.master) then
#endif
      GRIBIT2 imap,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iimap,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(iimap),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'imap grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 2 riv storage
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          if (gdriv(i,j).gt.1.E-10) then
            gdriv0(i,j)=log10(gdriv(i,j))
          else
            gdriv0(i,j)=99
          endif
        enddo
      enddo
      GRIBIT2 gdriv0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,igdriv,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(igdriv),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'gdriv grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 3 riv discharge
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          rflow0(i,j)=rflow(i,j)*rtime
          if (rflow0(i,j).gt.1.E-10) then
            rflow0(i,j)=log10(rflow0(i,j))
          else
            rflow0(i,j)=99
          endif
        enddo
      enddo
      GRIBIT2 rflow0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irflow,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irflow),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'rflow grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 4 total runoff
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          roff(i,j)=roff(i,j)*rtime
        enddo
      enddo
      GRIBIT2 roff,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iroff,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iroff),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'roff grib make failed'
      endif
#ifdef MP
      endif
#endif
c
#ifdef ISOTOPE
c
c 5 riv storage isotope 1
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          if (gdriv1(i,j).gt.1.E-10) then
            gdriv0(i,j)=gdriv1(i,j)/gdriv(i,j)
          else
            gdriv0(i,j)=0.
          endif
        enddo
      enddo
      GRIBIT2 gdriv0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,igdrivi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(igdrivi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'gdrivi1 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 6 riv discharge isotope 1
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          if (rflow(i,j)*rtime.gt.1.E-10) then
            rflow0(i,j)=rflow1(i,j)/rflow(i,j)
          else
            rflow0(i,j)=0.
          endif
        enddo
      enddo
      GRIBIT2 rflow0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irflowi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &           ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irflowi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'rflowi1 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 7 total runoff isotope 1
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          roff1(i,j)=roff1(i,j)*rtime
        enddo
      enddo
      GRIBIT2 roff1,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iroffi1,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iroffi1),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'roffi1 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 8 riv storage isotope 2
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          if (gdriv(i,j).gt.1.E-10) then
            gdriv0(i,j)=gdriv2(i,j)/gdriv(i,j)
          else
            gdriv0(i,j)=0.
          endif
        enddo
      enddo
      GRIBIT2 gdriv0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,igdrivi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(igdrivi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'gdrivi2 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 9 riv discharge isotope 2
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          if (rflow(i,j)*rtime.gt.1.E-10) then
            rflow0(i,j)=rflow2(i,j)/rflow(i,j)
          else
            rflow0(i,j)=0.
          endif
        enddo
      enddo
      GRIBIT2 rflow0,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irflowi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &           ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(irflowi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'rflowi2 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 10 total runoff isotope 2
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          roff2(i,j)=roff2(i,j)*rtime
        enddo
      enddo
      GRIBIT2 roff2,lbmr,0,io2_,jo2_,maxbit,90,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iroffi2,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(iroffi2),iens,
     &            rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
     &           ,gr,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,gr)
      else
        print *,'roffi2 grib make failed'
      endif
#ifdef MP
      endif
#endif
c
#endif
c
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'grib riv file created.'
#ifdef MP
      endif
#endif
c
#ifdef MP
      if(mype.eq.master) then
#endif
      do j=1,jo2_
        do i=1,io2_
          rflow(i,j)=0.e0
          roff(i,j)=0.e0
#ifdef ISOTOPE
          rflow1(i,j)=0.e0
          roff1(i,j)=0.e0
          rflow2(i,j)=0.e0
          roff2(i,j)=0.e0
#endif
        enddo
      enddo
      close(nn)
#ifdef MP
      endif
#endif
c
#endif
c
#endif !! for RIVER
      return
      end
