
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAP22P lcap22p
#define LCAPP lcapp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNT2P lnt2p
#define LNTP lntp
#define LNUV lnuv
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NCPUS ncpus
#define NPES npes
#define NROW nrow
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define RD rd
#define RV rv
#define RERTH rerth
#define TWOJ1 twoj1
      subroutine gloopa
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gloop  computes dynamic non-linear tendency terms
c of temp. div. ln(ps)
c computes predicted values of vorticity and moisture
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-13
c
c abstract:
c   program  starts with spectral  coefficients temp.
c   of vorticity, divergence, specific humidity, and
c   ln((psfc).  converts them to the gaussian grid at each
c   latitude and calls fidi,  for the northern and southern
c   hemispheres at the same time.  after return from fidi
c   sr.  completes calculation of tendencies of temp. div. and lnps.
c   specific humidity, and vorticity are predicted by sr. sigvor
c   all input/output  is via commons.
c
c program history log:
c   91-03-06  joseph sela
c   99-01-04  hann-ming henry juang	merge loopb and loopr in
c   99-05-14  hann-ming henry juang	add mpi
c
c usage:    call gloopa
c
c attributes:
c   language: fortran 77
c
#include <paramodel.h>
#include <comfibm.h>
#ifndef DYNAMIC_ALLOC
      parameter(ncpus=ncpus_)
#endif
c  sof include
#ifdef MP
#define JCAP1S jcap1p_
#define LNT2S lnt2p_
#define LATG2S latg2p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAPS lcapp_
#define LCAP22S lcap22p_
#define LONF22S lonf22p_
#else
#define JCAP1S jcap1_
#define LNT2S lnt2_
#define LATG2S latg2_
#define LEVSS levs_
#define LEVHS levh_
#define LCAPS lonf_
#define LCAP22S lonf22_
#define LONF22S lonf22_
#endif
#ifdef NISLQ
#include <comnislq.h>
#endif
c
c syn(1, 0*levs_+1, lan)  ze
c syn(1, 1*levs_+1, lan)  di
c syn(1, 2*levs_+1, lan)  te
c syn(1, 3*levs_+1, lan)  rq
c syn(1, 4*levs_+1, lan)  uln
c syn(1, 5*levs_+1, lan)  vln
c syn(1, 6*levs_+1, lan)  dpdphi
c syn(1, 6*levs_+2, lan)  dpdlam
c syn(1, 6*levs_+3, lan)  q
c
c dyn(1, 0*levs_+1, lan)  d(t)/d(phi)
c dyn(1, 1*levs_+1, lan)  d(rq)/d(phi)
c dyn(1, 2*levs_+1, lan)  d(t)/d(lam)
c dyn(1, 3*levs_+1, lan)  d(rq)/d(lam)
c dyn(1, 4*levs_+1, lan)  d(u)/d(lam)
c dyn(1, 5*levs_+1, lan)  d(v)/d(lam)
c dyn(1, 6*levs_+1, lan)  d(u)/d(phi)
c dyn(1, 7*levs_+1, lan)  d(v)/d(phi)
c
c anl(1, 0*levs_+1, lan)  z     dqdt
c anl(1, 0*levs_+2, lan)  uu    dudt
c anl(1, 1*levs_+2, lan)  vv    dvdt
c anl(1, 2*levs_+2, lan)  y     dtdt
c anl(1, 3*levs_+2, lan)  rt    drdt
c

#ifndef NISLQ
      parameter(lots    =5*levs_+levh_+4,lotst=2*levs_+1,
#else
      parameter(lots    =5*levs_+4,lotst=2*levs_+1,
#endif
     &          ksz     =1,
     &          ksd     =ksz+levs_,
     &          kst     =ksd+levs_,
#ifndef NISLQ
     &          ksr     =kst+levs_,
     &          ksu     =ksr+levh_, kstb=ksr+levh_,
#else
     &          ksu     =kst+levs_, kstb=kst+levs_,
#endif
     &          ksv     =ksu+levs_,
     &          kspphi  =ksv+levs_,
     &          ksplam  =kspphi+1,
     &          ksp     =ksplam+1,
     &          ksplap  =ksp+1) 
#ifndef NISLQ
      parameter(lotd    =6*levs_+2*levh_,
#else
      parameter(lotd    =6*levs_,
#endif
     &          kdtphi  =1,
#ifndef NISLQ
     &          kdrphi  =kdtphi+levs_,
     &          kdtlam  =kdrphi+levh_,
     &          kdrlam  =kdtlam+levs_,
     &          kdulam  =kdrlam+levh_,
#else
     &          kdtlam  =kdtphi+levs_,
     &          kdulam  =kdtlam+levs_,
#endif
     &          kdvlam  =kdulam+levs_,
     &          kduphi  =kdvlam+levs_,
     &          kdvphi  =kduphi+levs_)
#ifndef NISLQ
      parameter(lota    =3*levs_+levh_+1,lotat=2*levs_,
#else
      parameter(lota    =3*levs_+1,lotat=2*levs_,
#endif
     &          kap     =1,
     &          kau     =kap+1,
     &          kav     =kau+levs_,
     &          kat     =kav+levs_
#ifndef NISLQ
     &,          kar     =kat+levs_)
#else
     &                             )
#endif
c
#ifndef NISLQ
      parameter(lotss   =5*LEVSS+LEVHS+4,lotsts=2*LEVSS+1,
#else
      parameter(lotss   =5*LEVSS+4,lotsts=2*LEVSS+1,
#endif
     &          kszs    =1,
     &          ksds    =kszs+LEVSS,
     &          ksts    =ksds+LEVSS,
#ifndef NISLQ
     &          ksrs    =ksts+LEVSS,
     &          ksus    =ksrs+LEVHS, kstbs=ksrs+LEVHS,
#else
     &          ksus    =ksts+LEVSS, kstbs=ksts+LEVSS,
#endif
     &          ksvs    =ksus+LEVSS,
     &          kspphis =ksvs+LEVSS,
     &          ksplams =kspphis+1,
     &          ksps    =ksplams+1,
     &          ksplaps =ksps+1)
#ifndef NISLQ
      parameter(lotds   =6*LEVSS+2*LEVHS,
#else
      parameter(lotds   =6*LEVSS,
#endif
     &          kdtphis =1,
#ifndef NISLQ
     &          kdrphis =kdtphis+LEVSS,
     &          kdtlams =kdrphis+LEVHS,
     &          kdrlams =kdtlams+LEVSS,
     &          kdulams =kdrlams+LEVHS,
#else
     &          kdtlams =kdtphis+LEVSS,
     &          kdulams =kdtlams+LEVSS,
#endif
     &          kdvlams =kdulams+LEVSS,
     &          kduphis =kdvlams+LEVSS,
     &          kdvphis =kduphis+LEVSS )

#ifndef NISLQ
      parameter(lotas   =3*LEVSS+LEVHS+1,
#else
      parameter(lotas   =3*LEVSS+1,
#endif
     &          lotats=2*LEVSS,
     &          kaps    =1,
     &          kaus    =kaps+1,
     &          kavs    =kaus+LEVSS,
     &          kats    =kavs+LEVSS
#ifndef NISLQ
     &,          kars    =kats+LEVSS)
#else
     &                              )
#endif
c local common
      common/comdhc/spdlat(levs_,LATG2S)
c ...................................................
c local array
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
      dimension uo1(lnt22_,levs_), vo1(lnt22_,levs_), to1(lnt22_,levs_),
     &          qo1(lnt22_,levh_), wo1(lnt22_,levs_), ps1(lnt22_)
      dimension uo2(lnt22_,levs_), vo2(lnt22_,levs_), to2(lnt22_,levs_),
     &          qo2(lnt22_,levh_), wo2(lnt22_,levs_), ps2(lnt22_)
      dimension au1(lnt22_,levs_), av1(lnt22_,levs_), at1(lnt22_,levs_),
     &          aq1(lnt22_,levh_), ap1(lnt22_)
      dimension au2(lnt22_,levs_), av2(lnt22_,levs_), at2(lnt22_,levs_),
     &          aq2(lnt22_,levh_), ap2(lnt22_)
      dimension dv1(lnt22_,levs_), dv2(lnt22_,levs_)
#ifdef REV_FRC
      dimension ac1(lnt22_,levs_), ac2(lnt22_,levs_)
#endif
#ifdef RLX_FRC
      dimension tscl1(lnt22_,levs_), tscl2(lnt22_,levs_)
#endif
c
#ifdef DBG
#ifndef CLM_CWF
      dimension hadq(lnt22_,levs_)
#endif
      dimension hadt(lnt22_,levs_)
      dimension vadt(lnt22_,levs_), vadq(lnt22_,levh_)
      dimension adbc(lnt22_,levs_)
      dimension rlxt(lnt22_,levs_), rlxq(lnt22_,levh_)
#endif
      real*4 dum4(lnt22_)
      character*80 fno
      integer idateg(4), klev
      logical first
      data first/.true./
      save first,
     +     uo1, vo1, to1, qo1, wo1, ps1, uo2, vo2, to2, qo2, wo2, ps2,
#ifdef REV_FRC
     +     ac1, ac2,
#endif
#ifdef ARM
#ifdef RLX_FRC
     +     tscl1,tscl2,
#endif
#endif
     +     dv1, dv2, au1, av1, at1, aq1, ap1, au2, av2, at2, aq2, ap2
      save hours,houre,hour1,hour2

#else                   /* not scm */

c................................................

#ifdef PSPLIT
#ifdef MP
#define LONF22S lonf22p_
#else
#define LONF22S lonf22_
#endif
#include <compsplit.h>
#endif

#ifdef MP
#include <commpi.h>
#include <compspec.h>
c
      dimension
     1 syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_),
     2 dyf(lonf22_,lotds,latg2p_),dgr(lonf22p_,lotd,latg2p_),
     3 anf(lonf22_,lotas,latg2p_),gra(lonf22p_,lota,latg2p_)
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c ........................................................
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus
#endif
      dimension scs2(JCAP1S),
     2 syn(LCAP22S,lotss,NCPUSS),syntop(2,JCAP1S,lotsts),
     3 dyn(LCAP22S,lotds,NCPUSS),
     4 anl(LCAP22S,lotas,NCPUSS),anltop(2,JCAP1S,lotats),
     5 flp(2,JCAP1S,lotas,NCPUSS),flm(2,JCAP1S,lotas,NCPUSS)
c#undef NCPUSS

#ifdef HYBRID
      integer nvcn
      real xvcn
c      real ak5(levp1_),bk5(levp1_)
#endif

#ifdef PSPLIT
#ifdef DG

#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
      dimension tgmxl(NCPUSS),igmxl(NCPUSS),kgmxl(NCPUSS)
      dimension tgmnl(NCPUSS),igmnl(NCPUSS),kgmnl(NCPUSS)
#endif
#ifdef DG3
#include <comgda.h>
#endif
#endif
#endif                  /* end SCM directive */


#ifdef NISLQ
#ifdef HYBRID
       real, dimension(LONF22S,levs_+1)           ::  pdot
#else
       real, dimension(LONF22S,levs_+1)           ::  dot
#endif
       real, dimension(LONF22S,levh_,LATG2S)::tmp1,tmp2,tmp3,slq
#endif /* NISLQ end */


c
#ifdef PSPLIT 
#ifdef DG3
      dimension gda(nwgda,kdgda,NCPUSS)
#endif
#ifdef RAS
      logical ras
      parameter (ras=.true.)
      parameter (nsphys=1)
      parameter (cp=cp_, alhl=hvap_, grav=g_, rgas=rd_)
      dimension sig(levs_+1), prj(levs_+1), prh(levs_)
     *,         fpk(levs_),   hpk(levs_),   sgb(levs_)
     *,         ods(levs_), rasal(levs_),  prns(levs_/2)
     *,      rannum(200)
c
       call setras(levs_, si, sl, del, cp, rgas, deltim
     *,            nsphys, thour
     *,            sig, sgb, prh, prj, hpk, fpk, ods, prns
     *,            rasal, lm, krmin, krmax, nstrp
     *,            ncrnd, rannum, afac, ufac)
#endif
#ifdef RASV2
      logical ras
      parameter (ras=.true.)
      parameter (nsphys=1)
      parameter (cp=cp_, alhl=hvap_, grav=g_, rgas=rd_, hfus=hfus_)
      real sig(levs_+1), prj(levs_+1), sgb(levs_),   rasal(levs_)
     *,    rannum(200,5), pdd, dsfc(LONF2S,LATG2S)
      integer mct
c
c  pdd is the lowest pressure above which downdraft is allowed
c  mct is a number of cloud types
c
      parameter (pdd=600.0, mct=15)
      integer krmin, krmax, kfmax, ncrnd, kctop(mct+1)
      common /rassav/sig, sgb, prj, rasal, krmin, krmax, kfmax
     *,              kctop, ncrnd, dsfc
c
        frain = .5
        if(inistp.eq.1) frain = 1.
        dtp= deltim / frain
c
c  check fhour and thour
c
        call setrasv2(levs_,  si,  sl, cp_, rd_, dtp, nsphys, fhour
     *,             sig,   sgb, prj
     *,             rasal, krmin, krmax, kfmax, ncrnd, rannum
     *,             mct,   kctop, deltim, dsfc, LONF2S*LATG2S)
#endif
c
#ifndef O3C
       if(ntrac_.ge.2) call seto3p(idate,fhour,no3p,no3l)
#endif
c
#endif
#ifdef MP
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
      jstr=latstr(mype)
      jend=latstr(mype)+latlen(mype)-1
      lons2=lonlen(mype)*2
      lats2=latlen(mype)
#else
      llstr=0
      llens=jcap1_
      lons2=lonf2_
      lats2=latg2_
#endif
c
#ifdef SCM
c=========================================================================
c SCM : read dynamic tendency terms (advection)
c=========================================================================
c
      if (first) then
        hour1 = thour
        hour2 = hour1 + dtbdy
      endif
      fv = rv_/rd_ - 1.0
c
      timer = curtime
      if (timer.eq.hour1) then
       nng = 31
       call fnam('scmbdy',6,thour,fno,ncho)
       open(nng,file=fno(1:ncho),form='unformatted',status='old')
       read(nng) thourg, idateg, dtv
#ifdef REV_FRC
       read(nng) uo1, vo1, at1, ac1, aq1, wo1, dv1, ps1
       read(nng) uo2, vo2, at2, ac2, aq2, wo2, dv2, ps2
#else
#ifdef RLX_FRC
#ifndef TOGA
       read(nng) uo1, vo1, to1, qo1, at1, aq1, wo1, dv1, ps1, tscl1
       read(nng) uo2, vo2, to2, qo2, at2, aq2, wo2, dv2, ps2, tscl2
#else
       read(nng) uo1, vo1, to1, qo1, at1, aq1, wo1, dv1, ps1
       read(nng) uo2, vo2, to2, qo2, at2, aq2, wo2, dv2, ps2
#endif
#else
       read(nng) uo1, vo1, at1, aq1, wo1, dv1, ps1
       read(nng) uo2, vo2, at2, aq2, wo2, dv2, ps2
#endif
#endif
       close(nng)
       hours = hour1
       houre = hour2
       hour1 = hour2
       hour2 = hour1 + dtbdy
      endif
      ftim1 = (timer-hours)/dtbdy
      ftim2 = (houre-timer)/dtbdy
      print 985, kdt,deltim,curtime,thour,ftim1,ftim2
985   format('GLOOPA ',i5,5f8.2)
c
c ... set wind & pressure ...
c
      do l = 1, lnt2_
       do k = 1, levs_
         x(l,k) = ftim2*uo1(l,k) + ftim1*uo2(l,k)
         w(l,k) = ftim2*vo1(l,k) + ftim1*vo2(l,k)
         vvel(l,k) = ftim2*wo1(l,k) + ftim1*wo2(l,k)
#ifdef CLM_CWF
         wdiv(l,k) = ftim2*dv1(l,k) + ftim1*dv2(l,k)
#endif
       enddo
       z(l) = ftim2*ps1(l) + ftim1*ps2(l)
      enddo
c
c ... advection of T & q ...
c
      do l = 1, lnt2_
       do k = 1, levs_
         y(l,k) = 0.0
#ifdef CLM_CWF
         hadq(l,k) = 0.0
#endif
#ifdef DBG
         hadq(l,k) = 0.0
         hadt(l,k) = 0.0
         vadt(l,k) = 0.0
         rlxt(l,k) = 0.0
         adbc(l,k) = 0.0
#endif
       enddo
       do k = 1, levh_
         rt(l,k) = 0.0
#ifdef DBG
         vadq(l,k) = 0.0
         rlxq(l,k) = 0.0
#endif
       enddo
      enddo
c
#ifndef REV_FRC
#ifdef DBG
      call gfidiu(lnt2_,vvel,te,rq,q,del,sl,rdel2,y,rt,vadt,adbc)
#else
      call gfidiu(lnt2_,vvel,te,rq,q,del,sl,rdel2,y,rt)
#endif
#endif
c
      do l = 1, lnt2_
       do k = 1, levs_
#ifdef RLX_FRC
#ifndef TOGA
         advtscl = ftim2*tscl1(l,k) + ftim1*tscl2(l,k) ! relaxation time scale (s)
#else
         advtscl = 3600.*24.  ! from Iacobellis & Somerville (2000)
#endif
#endif

         advq    = -1.0*(ftim2*aq1(l,k) + ftim1*aq2(l,k))
#ifdef RLX_FRC
#ifdef CLM_CWF
         hadq(l,k) = advq
#endif
         qobs = ftim2*qo1(l,k) + ftim1*qo2(l,k)
         advq = rt(l,k) + advq + (qobs-rm(l,k))/advtscl
#else
         advq = rt(l,k) + advq
#endif
#ifdef DBG
         vadq(l,k) = rt(l,k)
#ifdef RLX_FRC
         rlxq(l,k) = (qobs-rm(l,k))/advtscl
#endif
#endif
#define QMIN qmin
         rt(l,k) = rm(l,k) + 2.0*deltim*advq
         qmin = QMIN
         rt(l,k) = max(rt(l,k),qmin)
       enddo
c
       if (ntotal_.gt.1) then
        do k = levs_+1, levh_
cld3v0         rt(l,k) = 0.0      ! no time change -> error
cld3v1         rt(l,k) = rm(l,k)  ! zero advection (horizontal & vertical)
          advq = rt(l,k)
#ifdef WSM2
          lev_qrs = levs_*(ntrac_+1) + 1
          if (k.ge.lev_qrs) advq = 0.0
#endif
          rt(l,k) = rm(l,k) + 2.0*deltim*advq ! only vertical advection
          qmin = QMIN
          rt(l,k) = max(rt(l,k),qmin)
        enddo
       endif
c
       do k = 1, levs_
#ifdef RLX_FRC
#ifndef TOGA
         advtscl = ftim2*tscl1(l,k) + ftim1*tscl2(l,k) ! relaxation time scale (s)
#else
         advtscl = 3600.*24.  ! from Iacobellis & Somerville (2000)
#endif
#endif
         adtv = -1.0*(ftim2*at1(l,k) + ftim1*at2(l,k))  ! advection of Tv
#ifdef REV_FRC
         adia = ftim2*ac1(l,k) + ftim1*ac2(l,k)
         advt = y(l,k) + adtv + adia
#else
#ifdef RLX_FRC
         tobs = ftim2*to1(l,k) + ftim1*to2(l,k)
         advt = y(l,k) + adtv + (tobs-tem(l,k))/advtscl
#else
         advt = y(l,k) + adtv
#endif
#endif
#ifdef DBG
         hadt(l,k) = adtv
#ifdef RLX_FRC
         rlxt(l,k) = (tobs-tem(l,k))/advtscl
#endif
#ifdef REV_FRC
         adbc(l,k) = adia
#endif
#endif
         y(l,k) = tem(l,k) + 2.0*deltim*advt
       enddo
      enddo
c
#ifdef DBG
      print 986, kdt,thour,hours,houre,ftim1,ftim2,dtbdy
986   format('BDY ',i5,6f8.3)
987   format(6e13.5)
#ifdef REV_FRC
      print 984,uo1(1,1),vo1(1,1),at1(1,1),ac1(1,1),aq1(1,1),
     +          wo1(1,1),exp(ps1(1))
      print 984,uo2(1,1),vo2(1,1),at2(1,1),ac2(1,1),aq2(1,1),
     +          wo2(1,1),exp(ps2(1))
984   format(7e13.5)
      print 987,x(1,1),w(1,1),y(1,1),rt(1,1),vvel(1,1),exp(z(1))
#else
      print 987,uo1(1,1),vo1(1,1),at1(1,1),aq1(1,1),wo1(1,1),exp(ps1(1))
      print 987,uo2(1,1),vo2(1,1),at2(1,1),aq2(1,1),wo2(1,1),exp(ps2(1))
      print 987,x(1,1),w(1,1),y(1,1),rt(1,1),vvel(1,1),exp(z(1))
#endif
c
      do k = 1, levs_
        write(95,988)
     +  curtime,kdt,k,x(1,k),w(1,k),y(1,k),rt(1,k),vvel(1,k),
     +  z(1),hadq(1,k),vadq(1,k),rlxq(1,k),
     +  hadt(1,k),vadt(1,k),adbc(1,k),rlxt(1,k)
988     format(f8.3,2i5,14e13.5)
      enddo
#endif
c
      first=.false.
c
#else
c=========================================================================
c GSM : calculate dynamic tendency terms
c=========================================================================
c
#ifdef MP
#ifndef NISLQ
      klev=3+ntotal_  ! ze,di,te,rq
#else
      klev=3          ! ze,di,te
#endif      
      call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,klev)
      call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define DIS dia
#define ZES zea
#define ULNS ulna
#define VLNS vlna
#else
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define DIS di
#define ZES ze
#define ULNS uln
#define VLNS vln
#endif
      call dellnp(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,
     1            llstr,llens,lwvdef)
      call deldfsp(QS,QLAPS,llstr,llens,lwvdef)
      call dzuvle(DIS,ZES,ULNS,VLNS,syntop(1,1,1),
     1            syntop(1,1,LEVSS+1),llstr,llens,lwvdef)

#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
c
      do n=1,2
        do j=1,JCAP1S
          do l=1,lotats
            anltop(n,j,l)=0.0
          enddo
        enddo
      enddo
#ifdef PSPLIT
#ifdef DG
      tgmx=-1.e20
      tgmn= 1.e20
#endif
#endif
c
#ifdef MP
#define LNT2X lln22p_
#define UUS uua
#define VVS vva
#define YS ya
#define RTS rta
#define ZS za
#else
#define LNT2X lnt2_
#define UUS uu
#define VVS vv
#define YS y
#define RTS rt
#define ZS z
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(uu,uua,vv,vva,y,ya,rt,rta),local(j,k)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(uu,uua,vv,vva,y,ya,rt,rta)
cmic$1 private(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
      do k=1,LEVSS
       do j=1,LNT2X
        UUS(j,k)=0.0
        VVS(j,k)=0.0
         YS(j,k)=0.0
       enddo
      enddo
#ifndef NISLQ
      do k=1,LEVHS
       do j=1,LNT2X
        RTS(j,k)=0.0
       enddo
      enddo
#endif
      do j=1,LNT2X
       ZS(j)=0. e 0
      enddo
#undef LNT2X
#undef UUS
#undef VVS
#undef YS
#undef RTS
#undef ZS
c
#ifndef MP
c compute latitude band limits
      last=mod(latg2_,ncpus)
      nggs=(latg2_-last)/ncpus
      if(last.ne.0)nggs=nggs+1
      inclat=ncpus
      lat1=1-ncpus
      lat2=0
      latdon=0
c
      do 10000 ngg=1,nggs
      if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
      lat1=lat1+ncpus
      lat2=lat2+inclat
#endif
c
#ifdef PSPLIT
#ifdef DG3

#ifdef MP
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif
      do lat=lat1,lat2
        lan=lat-latdon
        call getdia(LATX,nwgda*kdgda,gda(1,1,lan))
      enddo
#undef LATX
#endif

#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
      do k=1,NCPUSS
        tgmxl(k)=tgmx
        tgmnl(k)=tgmn
      enddo
#undef NCPUSS
#endif
#endif
c
#ifdef MP
      lat1=1
      lat2=latg2_
      latdon=0
#define ZES zea
#define TES tea
#else
#define ZES ze
#define TES te
#endif
c first lat loop
#ifdef ORIGIN_THREAD
c$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
c$&        colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
c$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
cmic$1 shared(colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
cmic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
c
      do 1000 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif
         
        call sums2i(ZES,syn(1,1,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,lotss)

        call sumtop(syn(1,kstbs,lan),syntop,qvv(1,lat),
     &              llstr,llensd,lwvdef,lotsts)
1000  continue
c
c compute merid. derivs. of temp. and moisture using qdd.
c
#ifdef ORIGIN_THREAD
c$doacross share(dyn,qdd,lat1,lat2,latdon,tea,te,
c$&        llstr,llens,lwvdef,lcapd,lcapdp,mype),
c$&        local(lat,lan,i,k,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(dyn,qdd)
cmic$1 shared(lat1,lat2,latdon,llens)
cmic$1 shared(tea,te,llstr,llens,lwvdef,lcapd,lcapdp,mype)
cmic$1 private(lat,lan,i,k,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k,llensd)
#endif
c
      do 1100 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif

#ifndef NISLQ
       klev=LEVSS+LEVHS  ! te,rq
#else
       klev=LEVSS        ! te
#endif
        call sums2i(TES,dyn(1,kdtphis,lan),qdd(1,lat),
     1              llstr,llensd,lwvdef,klev)
c   d(t)/d(phi)  d(rq)/d(phi) in s. hemi.
        do k=1,LEVSS
         do i=1,llensd*2
          dyn(i+LCAPS,kdtphis-1+k,lan)=-dyn(i+LCAPS,kdtphis-1+k,lan)
         enddo
        enddo
#ifndef NISLQ
        do k=1,LEVHS
         do i=1,llensd*2
          dyn(i+LCAPS,kdrphis-1+k,lan)=-dyn(i+LCAPS,kdrphis-1+k,lan)
         enddo
        enddo
#endif
1100  continue

#ifdef ORIGIN_THREAD
c$doacross share(dyn,rcs2,syn,lat1,lat2,latdon,
c$&              llstr,llens,lwvdef,lcapd,lcapdp,mype),
c$&              local(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(dyn,rcs2,llstr,llens,lwvdef)
cmic$1 shared(syn,lat1,lat2,latdon,lcapd,lcapdp,mype)
cmic$1 private(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
c
      do 1200 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif
c
c   calculate t rq u v zonal derivs. by multiplication with i*l
c
        do l=1,llensd
          j=lwvdef(llstr+l)
          scs2(l)=float(j)/rerth_ *rcs2(lat)
        enddo
c
c   calculate t rq u v zonal derivs. by multiplication with i*l
c   note scs2=rcs2*l/rerth_
c
c   d(t)/d(lam)  d(rq)/d(lam)  d(u)/d(lam)  d(v)/d(lam) .......
c


#ifndef NISLQ
        klev=3*LEVSS+LEVHS
#else
        klev=3*LEVSS
#endif
        do k=1, klev
          kdd=kdtlams-1+k
          kss=ksts-1+k
          do j=1,llensd
            dyn(       2*j-1,kdd,lan)=-syn(       2*j  ,kss,lan)*scs2(j)
            dyn(       2*j  ,kdd,lan)= syn(       2*j-1,kss,lan)*scs2(j)
            dyn(LCAPS+2*j-1,kdd,lan)=-syn(LCAPS+2*j  ,kss,lan)*scs2(j)
            dyn(LCAPS+2*j  ,kdd,lan)= syn(LCAPS+2*j-1,kss,lan)*scs2(j)
          enddo
        enddo
1200  continue
c
#undef TES
#undef ZES
c
#ifdef MP
      call mpnl2ny(syn,lcap22p_,latg2_,
     1             syf,lonf22_,latg2p_,lotss,1,lotss)
#ifndef NISLQ
           klev=4*levsp_+2*levhp_
#else
           klev=4*levsp_
#endif
      call mpnl2ny(dyn,lcap22p_,latg2_,
     1             dyf,lonf22_,latg2p_,lotds,1,klev)
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define SYNS syf
#define DYNS dyf
#else
#define SYNS syn
#define DYNS dyn
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(syf,syn,dyf,dyn,lat1,lat2,latdon,latdef,lcapd,lonfd),
c$&        local(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syf,syn,dyf,dyn)
cmic$& shared(lat1,lat2,latdon,latdef,lcapd,lonfd)
cmic$& private(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,k,j,lcapf,lonff)
#endif
       do 1300 lat=lat1,lat2
#ifdef REDUCE_GRID
         lcapf=lcapd(latdef(lat))
         lonff=lonfd(latdef(lat))
#else
         lcapf=jcap1_
         lonff=lonf_
#endif
         lan=lat-latdon
         call ftilonf(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,
     1                lcapf,lonff,latdef(lat),1)


#ifndef NISLQ
       klev=4*LEVSS+2*LEVHS
#else
       klev=4*LEVSS
#endif

         call ftilonf(DYNS(1,1,lan),DYNS(1,1,lan),2*klev,
     1                lcapf,lonff,latdef(lat),1)
         do k=1,LEVSS
           do j=1,lonff*2
             DYNS(j,kduphis-1+k,lan)=
     1       DYNS(j,kdvlams-1+k,lan)-SYNS(j,kszs-1+k,lan)
             DYNS(j,kdvphis-1+k,lan)=
     1      -DYNS(j,kdulams-1+k,lan)+SYNS(j,ksds-1+k,lan)
           enddo
         enddo
 1300  continue
c
#undef SYNS
#undef DYNS
c
#ifdef MP
c
#ifndef NISLQ
           klev=5+ntotal_
#else
           klev=5
#endif
       call mpnk2nx(syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,
     &              klev)
       call mpx2nx (syf,lonf22_,lotss,
#ifdef PSPLIT
     &              grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
#else
     &              grs,lonf22p_,lots,latg2p_,kspphis,kspphi,3)
#endif

#ifndef NISLQ
           klev=6+2*ntotal_
#else
           klev=6
#endif

       call mpnk2nx(dyf,lonf22_,lotds,
     &              dgr,lonf22p_,lotd,latg2p_,levsp_,levs_,1,1,
     &              klev) 
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define SYNS grs
#define DYNS dgr
#define ANLS gra
#define LATX lan
#else
#define SYNS syn
#define DYNS dyn
#define ANLS anl
#define LATX lat
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(dgr,dyn,rcs2,grs,syn,gra,anl,lat1,lat2,latdon,spdlat,
c$&              rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype),
c$&              local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(dyn,dgr,rcs2)
cmic$1 shared(syn,grs,anl,gra,lat1,lat2,latdon,spdlat)
cmic$1 shared(rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype)
cmic$1 private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
c
      do 2000 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(lan,mype)*2
#else
        lonsd2=lonfd(latdef(lat))*2
#endif
#else
        lonsd2=lons2
#endif
c
#ifdef MP
        if( lonsd2.gt.0 ) then
#endif
        do k=1,levs_
         do j=1,lonsd2
#ifndef HYBRID
          SYNS(j,kst-1+k,lan)=SYNS(j,kst-1+k,lan)-tov(k)
#endif
#ifdef PSPLIT
          ANLS(j,kau-1+k,lan)=dudtm(j,k,LATX)
          ANLS(j,kav-1+k,lan)=dvdtm(j,k,LATX)
          ANLS(j,kat-1+k,lan)=dtdtm(j,k,LATX)
#else
          ANLS(j,kau-1+k,lan)=0.
          ANLS(j,kav-1+k,lan)=0.
          ANLS(j,kat-1+k,lan)=0.
#endif
         enddo
        enddo
        do k=1,levh_
         do j=1,lonsd2
#ifdef NISLQ
           slq(j,k,lan)=slq_q2(j,k,LATX)
#else
#ifdef PSPLIT
          ANLS(j,kar-1+k,lan)=drdtm(j,k,LATX)
#else
          ANLS(j,kar-1+k,lan)=0.0
#endif
#endif
         enddo
        enddo
        do j=1,lonsd2
         ANLS(j,kap,lan)=0.0
        enddo
c
#ifdef HYBRID
c        do k=1,levp1_
c           ak5(k)=0.
c           bk5(k)=si(levp1_-k+1)
cc           bk5(k)=si(k)
c        enddo
        call gfidi_hyb(lonsd2,lan,
     $            ak5,bk5,
     1            SYNS(1,ksd,lan),SYNS(1,kst,lan),
     1            SYNS(1,ksz,lan),SYNS(1,ksu,lan),
     1            SYNS(1,ksv,lan),
#ifndef NISLQ
     1            SYNS(1,ksr,lan),
#else
     1            slq(1,1,lan),
#endif
     1            SYNS(1,kspphi,lan),SYNS(1,ksplam,lan),
     $            SYNS(1,ksp,lan),
     1            rbs2(LATX),
     1            spdlat(1,LATX),deltim,nvcn,xvcn,
     1            DYNS(1,kdtphi,lan),DYNS(1,kdtlam,lan),
#ifndef NISLQ
     1            DYNS(1,kdrphi,lan),DYNS(1,kdrlam,lan),
#else
     1            tmp1(1,1,lan),tmp2(1,1,lan),
#endif
     1            DYNS(1,kdulam,lan),DYNS(1,kdvlam,lan),
     1            DYNS(1,kduphi,lan),DYNS(1,kdvphi,lan),
     1            ANLS(1,kap,lan),ANLS(1,kat,lan),
#ifndef NISLQ
     1            ANLS(1,kar,lan),ANLS(1,kau,lan),
     1            ANLS(1,kav,lan))
#else
     1            tmp3(1,1,lan),ANLS(1,kau,lan),
     1            ANLS(1,kav,lan),pdot)
#endif


#else
        call gfidiu(lonsd2,
     1            SYNS(1,ksd,lan),SYNS(1,kst,lan),
     1            SYNS(1,ksz,lan),SYNS(1,ksu,lan),
     1            SYNS(1,ksv,lan),
#ifndef NISLQ
     1            SYNS(1,ksr,lan),
#else
     1            slq(1,1,lan),
#endif
     1            SYNS(1,kspphi,lan),SYNS(1,ksplam,lan),
     1            rbs2(LATX),del,rdel2,ci,p1,p2,h1,h2,tov,
     1            spdlat(1,LATX),
     1            DYNS(1,kdtphi,lan),DYNS(1,kdtlam,lan),
#ifndef NISLQ
     1            DYNS(1,kdrphi,lan),DYNS(1,kdrlam,lan),
#else
     1            tmp1(1,1,lan), tmp2(1,1,lan), 
#endif
     1            DYNS(1,kdulam,lan),DYNS(1,kdvlam,lan),
     1            DYNS(1,kduphi,lan),DYNS(1,kdvphi,lan),
     1            ANLS(1,kap,lan),ANLS(1,kat,lan),
#ifndef NISLQ
     1            ANLS(1,kar,lan),ANLS(1,kau,lan),
     1            ANLS(1,kav,lan))
#else
     1            tmp3(1,1,lan),ANLS(1,kau,lan),
     1            ANLS(1,kav,lan),dot)
#endif


#endif

#ifdef NISLQ
c      dyn_grid for nislq
c
c      Psfc
c
       do i = 1,lonsd2
         slq_psfc2(i,LATX)=exp(SYNS(i,ksp,lan))
       enddo
c
c      U, V
c
       do k = 1,levs_
         do i=1,lonsd2
           slq_u2(i,levs_+1-k,LATX)=SYNS(i,ksu+k-1,lan)
           slq_v2(i,levs_+1-k,LATX)=SYNS(i,ksv+k-1,lan)
         enddo
       enddo
c
c      Pdot
c
       do k = 1,levs_+1
         do i = 1,lonsd2
#ifdef HYBRID
c
c          pdot:top to bottom
c
           slq_w2(i,k,LATX)=pdot(i,k)
#else
c
c          dot:bottom to top
c
           slq_w2(i,levs_+2-k,LATX)=-dot(i,k)*exp(SYNS(i,ksp,lan))
#endif
         enddo
       enddo

#endif /* NISLQ end */





#ifdef PSPLIT
c
#ifndef HYBRID
        do k=1,levs_
          do j=1,lonsd2
            SYNS(j,kst-1+k,lan)=SYNS(j,kst-1+k,lan)+tov(k)
          enddo
        enddo
#endif
        if(.not.stepone) then
          call gbphys(lonsd2,
     x              SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),
     x              SYNS(1,ksu,lan),SYNS(1,ksv,lan),
     x              SYNS(1,ksp,lan),SYNS(1,kst,lan),
     x              SYNS(1,ksr,lan),SYNS(1,ksd,lan),
     x              SYNS(1,ksplap,lan),
     x              dtdtm(1,1,LATX),drdtm(1,1,LATX),
     x              dudtm(1,1,LATX),dvdtm(1,1,LATX),
#ifdef DG
     x              tgmxl(lan),igmxl(lan),kgmxl(lan),
     x              tgmnl(lan),igmnl(lan),kgmnl(lan),
#endif
#ifdef DG3
     x              gda(1,1,lan),
#endif
#ifdef RAS
     x              ras,levs_,cp,alhl,grav,rgas,
     x              sig, sgb, prh, prj, hpk, fpk, ods, prns,
     x              rasal, lm, krmin, krmax, nstrp,
     x              ncrnd, rannum, afac, ufac,
#endif
#ifdef RASV2
     &              ras,rgas, cp, grav, alhl,
     &              sig,prj,sgb,rasal, rannum, dsfc(1,LATX),
     &              pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,
#endif
#ifndef VIC
     x              LATX,1.0)
#else
     x              LATX,1.0,idate)
#endif
        endif
c
#ifdef DG
        if(tgmxl(lan).gt.tgmx) then
          tgmx=tgmxl(lan)
          igmx=igmxl(lan)
          kgmx=kgmxl(lan)
          jgmx=lat
        else if(tgmnl(lan).lt.tgmn) then
          tgmn=tgmnl(lan)
          igmn=igmnl(lan)
          kgmn=kgmnl(lan)
          jgmn=lat
        endif
#endif
#endif
#ifdef MP
        endif
#endif
c
2000  continue
#ifdef PSPLIT
#ifdef DG3
#ifdef MP
#define LATX lan
#else
#define LATX lat
#endif
      do lat=lat1,lat2
        lan=lat-latdon
        call putdia(LATX,nwgda*kdgda,gda(1,1,lan))
      enddo
#endif
#endif
#undef SYNS
#undef DYNS
#undef ANLS
#undef LATX
c
#ifdef MP
#ifndef NISLQ
       klev=3+ntotal_
#else
       klev=3
#endif
      call mpnx2nk(gra,lonf22p_,lota,
     &             anf,lonf22_ ,lotas,latg2p_,levs_,levsp_,kau,kaus,
     &             klev) 
      call mpnx2x (gra,lonf22p_,lota,
     &             anf,lonf22_ ,lotas,latg2p_,kap,kaps,1)
      if( levlen(mype).lt.levsp_ ) then
        do lan=1,lats2
         do i=1,lonf22_
          anf(i,kaus+levsp_-1,lan)=0.0
          anf(i,kavs+levsp_-1,lan)=0.0
          anf(i,kats+levsp_-1,lan)=0.0
         enddo
        enddo
#ifndef NISLQ
        do lan=1,lats2
         do i=1,lonf22_
           do kc = 1,ntotal_
             anf(i,kars+levsp_*kc-1,lan)=0.0
           enddo
         enddo
        enddo
#endif
      endif
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define ANLS anf
#else
#define ANLS anl
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd),
c$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd)
cmic$1 private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
c
      do 2200 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
        lcapf=lcapd(latdef(lat))
        lonff=lonfd(latdef(lat))
#else
        lcapf=jcap1_
        lonff=lonf_
#endif
        call ftilonf(ANLS(1,1,lan),ANLS(1,1,lan),2*lotas,
     1               lcapf,lonff,latdef(lat),-1)
2200  continue
#undef ANLS
c
#ifdef MP
      call mpny2nl(anf,lonf22_ ,latg2p_,
     1             anl,lcap22p_,latg2_ ,lotas,kaps,lotas)
c
      lat1=1
      lat2=latg2_
      latdon=0
#define ZS za
#else
#define ZS z
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype),
c$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype)
cmic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
c
      do 2500 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif
c
        call flpflm(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),
     1            llensd,lotas)
2500  continue
c 
c   no multi-threads should be given in following loop
c   otherwise results will be non-produceable
c
      do 3000 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif
c
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
        call fl2i(flp(1,1,1,lan),flm(1,1,1,lan),ZS,qww(1,lat),
     1            llstr,llensd,lwvdef,lotas)
#endif
#ifdef DEFAULT
        call fl2ip(flp(1,1,1,lan),flm(1,1,1,lan),ZS,qww(1,lat),
     1             llstr,llensd,lwvdef,lotas)
#endif
c
        call uvsums(flp(1,1,kaus,lan),flm(1,1,kaus,lan),
     1              flp(1,1,kavs,lan),flm(1,1,kavs,lan),
     2              anltop(1,1,1),anltop(1,1,LEVSS+1),
     3              qvv(1,lat),wgt(lat),
     4              llstr,llensd,lwvdef,LEVSS)
c
3000  continue
#undef ZS
c
#ifndef MP
      latdon=latdon+(lat2-lat1+1)
10000 continue
c
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(spdmax,spdlat,lats2),local(k,lat)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(spdmax,spdlat,lats2)
cmic$1 private(k,lat)
#endif
#ifdef OPENMP
!$omp parallel do private(k,lat)
#endif
      do k=1,levs_
        spdmax(k) = 0.0
        do lat=1,lats2
          spdmax(k)=max(spdmax(k),spdlat(k,lat))
        enddo
        spdmax(k)=sqrt(spdmax(k))
      enddo
#ifdef MP
      call mpgetspd(spdmax)
      if( mype.eq.master ) then
#endif
#ifndef NOPRINT
      print 100,(spdmax(k),k=1,levs_)
100   format(' global checked speed maxima for all layers ',
     &     :/' spdmx(01:10)=',10f5.0,:/' spdmx(11:20)=',10f5.0,
     &     :/' spdmx(21:30)=',10f5.0,:/' spdmx(31:40)=',10f5.0,
     &     :/' spdmx(41:50)=',10f5.0,:/' spdmx(51:60)=',10f5.0,
     &     :/' spdmx(61:70)=',10f5.0,:/' spdmx(71:80)=',10f5.0,
     &     :/' spdmx(81:90)=',10f5.0,:/' spdmx(91:00)=',10f5.0)
#endif
#ifdef MP
      endif
#endif
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      do k=1,levs_
        if(spdmax(k).eq.0.) then
          print *,'run failure.  spdmax=0.'
          call abort
        endif
      enddo
#ifdef MP
      endif
#endif
c
c     input : w=d(u)/d(t) x=d(v)/d(t)
c     output: uln=d(di)/d(t) vln=d(ze)/d(t)
c
#ifdef MP
#define UUS uua
#define VVS vva
#define XS xa
#define WS wa
#else
#define UUS uu
#define VVS vv
#define XS x
#define WS w
#endif
      call uvtodz(UUS,VVS,XS,WS,anltop(1,1,1),
     1            anltop(1,1,LEVSS+1),llstr,llens,lwvdef)
      
   
#undef UUS
#undef VVS
#undef XS
#undef WS
c
#ifdef MP
#ifndef NISLQ
      klev=3+ntotal_  ! y,rt,x,w
#else
      klev=3          ! y,x,w
#endif
      call mpnk2nn(ya,lln22p_,levsp_, y,lnt22p_,levs_,klev)
      call mpn2nn(za,lln22p_,z,lnt22p_,1)
#endif
c
c  subtract off linear dependence on divergence
c
      do k=1,levs_
        do j=1,levs_
          do i=1,LNT2S
#ifndef HYBRID
            y(i,k)=y(i,k)-bm(k,j)*di(i,j)
#else
            y(i,k)=y(i,k)+bm(k,j)*di(i,j)
#endif
          enddo
        enddo
      enddo
#ifdef HYBRID
      do j=1,levs_
        do i=1,LNT2S
          z(i)=z(i)+sv(j)*di(i,j)
        enddo
      enddo
#endif
c
c add topog. contriba into x.
c integrate vorticity amd moisture in time
c
      do k=1,levs_
       do i=1,LNT2S
         x(i,k)=x(i,k)+gz(i)
         w(i,k)=zem(i,k)+2.*deltim*w(i,k)
       enddo
      enddo

#ifndef NISLQ
      do k=1,levh_
       do i=1,LNT2S
        rt(i,k)= rm(i,k)+2.*deltim* rt(i,k)
       enddo
      enddo
#endif /* NISLQ end */

#ifdef MP
      if( mype.eq.master ) then
#endif
      do k=1,levs_
       w(1,k)=0.e0
       w(2,k)=0.e0
      enddo
#ifdef MP
      endif
#endif
#ifdef PSPLIT
#ifdef DG
#ifdef MP
      if(mype.eq.master) then
#endif
        print *,'(" gloopb t range ",2(4x,f6.1," @i,k,lat ",3i4))',
     &        tgmx,igmx,kgmx,jgmx,tgmn,igmn,kgmn,jgmn
#ifdef MP
      endif
#endif
#endif
#endif

#endif
c
      return
      end
