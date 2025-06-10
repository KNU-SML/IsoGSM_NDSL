#define BGF bgf
#define BORDER border
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
#ifndef RSM
      subroutine gradpre(idate,dorad)
#else
      subroutine gradpre(dorad)
#endif
#include <paramodel.h>
#ifndef RSM
#include <comfgrid.h>
#include <comfphys.h>
#include <comfver.h>
#include <comio.h>
#include <comgrad.h>
#include <radiag.h>
#ifdef MP
#include <commpi.h>
      dimension cvf  (lonf2_,latg2_)
     &         ,cvtf (lonf2_,latg2_)
     &         ,cvbf (lonf2_,latg2_)
     &         ,cvrf (lonf2_,latg2_)
     &         ,cvtrf(lonf2_,latg2_)
     &         ,cvbrf(lonf2_,latg2_)
#define LATG2S latg2p_
#define LONF2S lonf2p_
#else
#define LATG2S latg2_
#define LONF2S lonf2_
#endif
#else
#include <rscomf.rerun.h>
#include <rscomltb.h>
#include <rscommap.h>
#include <rscomgrad.h>
#include <rdparm.h>
#ifdef MP
#include <comrmpi.h>
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif
c
      parameter (nozon=48)
      parameter (loz=17)
      common /sbuv/ psnasa(loz),o3nasa(37,loz)
      dimension pstr(loz)
c
#ifndef RSM
      dimension idate(4)
#endif
      logical dorad
      integer season
      common   /diucon/
     2   season,             lseason,               fcstda,
     3   jtyme(5),           lftype, daz(12),    jdnmc, ljdnmc,
     4   fjdnmc,             tslag,              rlag,
     5   timin,              tpi,                hpi,
     6   year,               day,                dhr,
     7   ixxxx,            lixxxx
c
      common /tables/  sc
c
      parameter(lots=5*levs_+3)
      parameter(lworkr=(20/lots*20+lots/20*lots)/(20/lots+lots/20))
#ifndef RSM
      dimension workr(lonf2_,lworkr)
      dimension iworkr(lonf2_)
#else
      dimension workr(igrd12_,lworkr)
      dimension iworkr(igrd12_)
      dimension latdef(jgrd12_)
c
      do k=1,jgrd12_
        latdef(k)=k
      enddo
#endif
c
      dthr = deltim / 3600.0
      hdthr = 0.5 * dthr
c
      itimsw = 0
      itimlw = 0
      dtsmod = amod(solhr,dtswav)
      if(inistp.ne.0 .or. dtsmod.lt.hdthr .or. dtsmod.ge.dtswav-hdthr)
     1   itimsw = 1
      dtlmod = amod(solhr,dtlwav)
      if(inistp.ne.0 .or. dtlmod.lt.hdthr .or. dtlmod.ge.dtlwav-hdthr)
     1   itimlw = 1
      dorad=.false.
      if(itimsw.eq.0 .and. itimlw.eq.0) return
      if(.not.runrad .and. inistp.ne.0) return
      dorad=.true.
c
cdg3     compute num sec between calls to cld code...for diagnostics
c
#ifdef DG
      dtacc = min(dtswav,dtlwav)*3600.
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print 1001, jcap_, levs_
 1001 format (1h0,'gfdl reduced rad',i2,i2,'g,e typ, feb 20 1986')
#ifdef MP
      endif
#endif
#endif
c
c  paerf when SWRMDC =.false.
c  alvsf,alnsf,alvwf,alnwf from the albaer.snl will be discarded
c   when NEWALB = .true.
c
c....
      if(inistp.ne.0) then
#ifndef RSM
        call insur2(nalaer,alvsf,alnsf,alvwf,alnwf,
#ifndef SW99UPDATE
     1             facsf,facwf,paerf)
#else
     1             facsf,facwf)
#endif
#else
        call rinsur2(LONF2S,LATG2S,nalaer,alvsf,alnsf,alvwf,alnwf,
#ifndef SW99UPDATE
     1             facsf,facwf,paerf)
#else
     1             facsf,facwf)
#endif
#endif
      endif
c
#ifndef RSM
      do j=1,LATG2S
        do i=1,LONF2S
          slmskr(i,j) =  slmsk(i,j)
        enddo
      enddo
#endif
c
c   astronomy calculations-once for each new radiation step
c 
c      get 4 digit year for julian day computation
c
      iyr = idate(4)
      imon = idate(2)
      iday = idate(3)
      iztim = idate(1)
      if(iyr.lt.100) then
        print *,'iy .lt.100 in gradpre'
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
      kyear = iyr
#ifndef RSM
      chour=fhour+shour/3600.
#else
      chour=fhour
      if(inistp.eq.0) chour=thour
#endif
c
c following broadcasts are to solve the bit reproducibility issue
c
      call compjd(kyear,imon,iday,iztim,0,jdnmc,fjdnmc)
      call fcstim(chour,imon,iday,iztim,jdnmc,fjdnmc,rlag,year,
     1            rsin1,rcos1,rcos2,jd,fjd)
      if(itimsw.eq.1) then
        call solar(jd,fjd,r1,dlt,alf,slag,sdec,cdec)
        solc=sc/(r1*r1)
        call coszmn(LONF2S,LATG2S,
#ifndef RSM
     1              dtswav,solhr,sinlab,coslab,sdec,cdec,slag,
     2              xlon,coszer,.true.,coszdg)
#else
     1              dtswav,solhr,sinlar,coslar,sdec,cdec,slag,
     2              rlon,coszer,.true.,coszdg)
#endif
      endif
c
      call cdate(jd,fjd,munth,im,id,iyear,ihr,xmin)
      call prtime(id,munth,iyear,ihr,xmin,jd,fjd,dlt,alf,r1,slag,solc)
c
c   nasa o3   calculations-once for each new radiation step      
c   get new climo from nasa 12-month sbuv data                  
c   if jo3=0 use old gfdl climo, jo3=1 use new nasa climo ozone
c
#ifndef O3C
      if(thour.eq.0.or.ntrac_.le.1) then
#endif
        jo3=0
        call o3sbuv(fhour,idate,nozon,o3nasa,pstr,jerr)
        if (jerr.le.0) then
c
c  if nasa data file was available (jerr=0), not (jerr=1)
c
          jo3=1
c
c   get nasa pressure in cb (flip vertical coordinate)
c
          do n=1,loz
            psnasa(n) = pstr(loz+1-n)*1.e-1
          enddo
        end if
#ifdef DBG
        if(mype.eq.master ) then
          print 167
          if (jo3.eq.0) print 166
          if (jo3.eq.1) print 168
          print 167
  166 format('    using gfdl zonal seasonal ozone climo ')
  167 format('       ------        ')
  168 format('    using nasa zonal monthly ozone climo ')
        endif
#endif

#ifndef O3C
      endif
#endif

#ifndef MP
#ifdef DBG
      print *,'debug print from gradpre'
      call maxmin(tsea  ,lonf2_,latg2_,1,1,'tsea')
      call maxmin(sheleg,lonf2_,latg2_,1,1,'sheleg')
      call maxmin(albedo,lonf2_,latg2_,1,1,'albedo')
      call maxmin(slmsk ,lonf2_,latg2_,1,1,'slmsk')
#endif
#endif
c
c  albedo and aerosol data interpolations
c
      if(itimsw.eq.1) then
#ifndef RSM
        do j=1,LATG2S
          do i=1,LONF2S
            coszen(i,j)=coszer(i,j)
          enddo
        enddo
#endif
#ifdef MP
      latstart=latstr(mype)
#else
      latstart=1
#endif
#ifdef ALB4
        call albaer(LONF2S,LATG2S,
#ifndef RSM
     1              im,slmsk,sheleg,zorl,coszen,tsea,hprime,jsno,
#else
     1              im,slmsk,sheleg,zorl,coszer,tsea,hprime,jsno,
#endif
     2              albedo(1,1,1),albedo(1,1,3),albedo(1,1,2),
     3              albedo(1,1,4),
#ifndef SW99UPDATE
#ifndef RSM
     3              facalf(1,1,1),facalf(1,1,2),paerf,
     4              xlat,alvbr,alnbr,alvdr,alndr,paerr,latdef(latstart))
#else
     3              facalf(1,1,1),facalf(1,1,2),paerf,
     4              rlat,alvbr,alnbr,alvdr,alndr,paerr,latdef(latstart))
#endif
#else
#ifndef RSM
     3              facalf(1,1,1),facalf(1,1,2),
     4              xlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)
#else
     3              facsf,facwf,
     4              rlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)
#endif
     5,             iswsrc(1),nfaer,kprfg,idxcg,cmixg,denng)
#endif
#else
        call albaer(LONF2S,LATG2S,
#ifndef RSM
     1              im,slmsk,sheleg,zorl,coszen,tsea,hprime,jsno,
#else
     1              im,slmsk,sheleg,zorl,coszer,tsea,hprime,jsno,
#endif
#ifndef SW99UPDATE
     2              alvsf,alnsf,alvwf,alnwf,facsf,facwf,paerf,
#ifndef RSM
     3              xlat,alvbr,alnbr,alvdr,alndr,paerr,latdef(latstart))
#else
     3              rlat,alvbr,alnbr,alvdr,alndr,paerr,latdef(latstart))
#endif
#else
     2              alvsf,alnsf,alvwf,alnwf,facsf,facwf,
#ifndef RSM
     3              xlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)
#else
     3              rlat,alvbr,alnbr,alvdr,alndr,latdef(latstart)
#endif
     4,             iswsrc(1),nfaer,kprfg,idxcg,cmixg,denng)
#endif
#endif
      endif
c
      tsmin=tsea(1,1)
      tsmax=tsea(1,1)
      shmin=sheleg(1,1)
      shmax=sheleg(1,1)
c
#ifdef RAD_SMOOTH_CLOUD
#ifdef MP
      call mpgp2f(cv ,lonf2p_,latg2p_,cvf ,lonf2_,latg2_,1)
      call mpgp2f(cvt,lonf2p_,latg2p_,cvtf,lonf2_,latg2_,1)
      call mpgp2f(cvb,lonf2p_,latg2p_,cvbf,lonf2_,latg2_,1)
      if( mype.eq.master ) then
#define CVX cvf
#define CVTX cvtf
#define CVBX cvbf
#define CVRX cvrf
#define CVTRX cvtrf
#define CVBRX cvbrf
#else
#define CVX cv
#define CVTX cvt
#define CVBX cvb
#define CVRX cvr
#define CVTRX cvtr
#define CVBRX cvbr
#endif
      call cvintf(CVX,CVTX,CVBX,lonf2_,latg2_,latg2_,
     1         CVRX,CVTRX,CVBRX,lonf2_,latg2_,latg2_,
     3         workr(1,1),workr(1,5),workr(1,9),workr(1,13),
     4         workr(1,17),iworkr,1,1,1)
#undef CVX
#undef CVTX
#undef CVBX
#undef CVRX
#undef CVTRX
#undef CVBRX
#ifdef MP
      endif
      call mpgf2p(cvrf ,lonf2_,latg2_,cvr ,lonf2p_,latg2p_,1)
      call mpgf2p(cvtrf,lonf2_,latg2_,cvtr,lonf2p_,latg2p_,1)
      call mpgf2p(cvbrf,lonf2_,latg2_,cvbr,lonf2p_,latg2p_,1)
#endif
#else
      do j=1,LATG2S
        do i=1,LONF2S
          cvr  (i,j)=cv   (i,j)
          cvtr (i,j)=cvt  (i,j)
          cvbr (i,j)=cvb  (i,j)
        enddo
      enddo
#endif
c
      return
      end
