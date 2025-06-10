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
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LNT2 lnt2
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NPES npes
#define PI pi
#include <define.h>
      subroutine nfdr_wriflx(nn,iyr,imo,ida,ihr,ifhr,thour)
#include <paramodel.h>
#ifndef RSM
#include <comio.h>
#endif

#ifdef RSM
#include <rscommap.h>
#endif
      parameter(it2pm=100,irh2pm=101,icld2pm=102,iw2pm=103,
     $          ipdur=104,irhmax=105,irhmin=106,ippt24=96,
     $          itmx24=97,itmn24=98)
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
#ifndef RSM
#define ILEN lonf_
#define JLEN latg_
#ifdef MP
#define MPGP2F mpgp2f
#define MPABORT mpabort
#include <commpi.h>
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else     /* MP */
#define MPABORT abort
#define LONF2S lonf2_
#define LATG2S latg2_
#endif      /* MP */
#else     /* RSM */
#define ILEN igrd1_
#define JLEN jgrd1_
#ifdef MP
#define MPGP2F rmpgp2f
#define MPABORT rmpabort
#include <comrmpi.h>
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else     /* MP */
#define MPABORT abort
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif      /* MP */
#endif      /* RSM */
c
#ifdef RSM
#include <rscomloc.h>
#endif
#include <comnfdr.h>
c
c  full arrays
c
      parameter(len=ILEN*JLEN)
      logical   lbm(len)
      character g(200+len*(32+1)/8)
c
      parameter(rad=pi_/180.)
c
      parameter(lens=ILEN2S*JLEN2S)
c
c  full working arrays
c
      dimension work(len)
c
      dimension ids(255)
      dimension iens(5)
c
      character*128 fno
c
      if(nn.le.0) return
c
      ithr=thour
#ifdef MP
      if(mype.eq.master) then
#endif
      call fnam('nfdr.frcng',10,thour,fno,ncho)
      open(unit=nn,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file ',fno(1:ncho)
      call abort
  901 continue
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
      rewind nn
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
      if(ithr.eq.0) ithr=1
c
#ifdef RSM
      proj=rproj
      delx=rdelx
      dely=rdely
      nproj=rproj
      rrlat1=rlat1*rad
      rrlat2=rlat2*rad
      rrlon1=rlon1*rad
      rrlon2=rlon2*rad
      if( nproj.eq.0 ) then   ! mercater
        idrt=1
        ortru=rtruth
        proj=0.0
      elseif( abs(nproj).eq.1 ) then  ! polar projection
        idrt=5
        ortru=rorient
        proj=rproj
      else
        print *,' error in nfdr_wriflx for projection type.'
      endif
#else
      idrt=4
      rrlat1=0.
      rrlon1=0.
      rrlat2=0.
      rrlon2=0.
      delx=0.
      dely=0.
      ortru=0.
      proj=0.
#endif
c
      maxbit=16
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
c
c  grib output for nfdr 
c
c 1.  2m temperature at 2pm
c
#ifdef MP
      call MPGP2F(t2m2pm,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=t2m2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'t2m2pm')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,it2pm,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(it2pm),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'t2m2pm grib created    '
#endif
      else
        print *,'t2m2pm grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 2. 2m relative humidity at 2pm
c
#ifdef MP
      call MPGP2F(rh2m2pm,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=rh2m2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rh2m2pm')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irh2pm,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(irh2pm),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'rh2m2pm grib created    '
#endif
      else
        print *,'rh2m2pm grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 3.  total cloudiness at 2pm
c
#ifdef MP
      call MPGP2F(cloud2pm,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=cloud2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'cloud2pm')
#endif
      il1k=0 
      il2k=0 
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,icld2pm,icolmn,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(icld2pm),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'cloud2pm grib created    '
#endif
      else
        print *,'cloud2pm grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 4.  10m wind speed at 2pm
c
#ifdef MP
      call MPGP2F(wind2pm,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=wind2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'wind2pm')
#endif
      il1k=0
      il2k=10
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,iw2pm,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(iw2pm),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'wind2pm grib created    '
#endif
      else
        print *,'wind2pm grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 5.  24 hour accumulated precipitation
c
#ifdef MP
      call MPGP2F(ppt24,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
        work(n)=ppt24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'ppt24')
#endif
      il1k=0
      il2k=0
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ippt24,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(ippt24),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'accumulated precip  grib created    '
#endif
      else
        print *,'accumulated precip  grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 6.  precipitation duration
c
#ifdef MP
      call MPGP2F(prcpdur,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=prcpdur(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'prcpdur')
#endif
      il1k=0
      il2k=0
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,ipdur,isfc,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,ina,inm,icen2,ids(ipdur),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'prcpdur grib created    '
#endif
      else
        print *,'prcpdur grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 7.  24 hour max 2m temperature
c
#ifdef MP
      call MPGP2F(tmax24,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=tmax24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'tmax24')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itmx24,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itmx24),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'24 hour tmax grib created    '
#endif
      else
        print *,'24 hour tmax make failed'
      endif
#ifdef MP
      endif
#endif
c
c 8.  24 hour min 2m temperature
c
#ifdef MP
      call MPGP2F(tmin24,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=tmin24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rhmin')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,itmn24,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(itmn24),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'24 hour min temperature  grib created    '
#endif
      else
        print *,'24 hour min temperature  grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 9.  max 2m relative humidity
c
#ifdef MP
      call MPGP2F(rhmax,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=rhmax(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rhmax')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irhmax,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(irhmax),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'rhmax grib created    '
#endif
      else
        print *,'rhmax grib make failed'
      endif
#ifdef MP
      endif
#endif
c
c 10.  min 2m relative humidity
c
#ifdef MP
      call MPGP2F(rhmin,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
         work(n)=rhmin(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(work,1)
#else
      call rowsep(work)
#endif
#ifdef DBG
      call maxmin(work,ILEN*JLEN,1,1,1,'rhmin')
#endif
      il1k=0
      il2k=2
      call gribit( work,lbm,idrt,ILEN,JLEN,maxbit,colat,
     &            ilpds,iptv,icen,igen,
     &            ibm0,irhmin,ielev,il1k,il2k,iyr,imo,ida,ihr,
     &            ifhour,ithr,0,inst,ina,inm,icen2,ids(irhmin),iens,
     &            rrlat1,rrlon1,rrlat2,rrlon2,delx,dely,ortru,proj
     &           ,g,lg,ierr)
      if(ierr.eq.0) then
        call wryte(nn,lg,g)
#ifdef DBG
        print *,'rhmin grib created    '
#endif
      else
        print *,'rhmin grib make failed'
      endif
#ifdef MP
      endif
#endif
c
#ifdef MP
      if(mype.eq.master) then
#endif
#ifndef NOPRINT
      print *,'grib nfdr forcing file created.'
#endif
      close(nn)
#ifdef MP
      endif
#endif
c
      return
      end
