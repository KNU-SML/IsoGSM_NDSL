#define IO2 io2
#define JO2 jo2
#define BGF bgf
#define BORDER border
#define CAL cal
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
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
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define MLVARK mlvark
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NP np
#define NQ nq
#define NSOIL nsoil
#define NT nt
#define NTOTAL ntotal
#define NTRAC ntrac
#define NU nu
#define NV nv
#define QMIN qmin
#define RD rd
#define RLPNT rlpnt
#define RLTSTP rltstp
#define RMLVARK rmlvark
#define RSLVARK rslvark
#define RV rv
#define SBC sbc
#define SLVARK slvark
#include <define.h>
#ifndef RSM
      subroutine gbphys(lons2,
#ifdef SCM
     x                  uugrs,vvgrs,plgr,tvgrs,rqgrs,xxgrs,
#else
     x                  pllamgr,plphigr,uugrs,vvgrs,plgr,
     x                  tvgrs,rqgrs,xxgrs,
     x                  pslap,
#endif
     x                  ggt0,ggq0,ggu0,ggv0,
#ifdef DG
     x                  tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,
#endif
#ifdef DG3
     x                  gda,
#endif
#else
      subroutine gbphys(
#ifdef SCM
     x                  ugrs,vgrs,pgr,tgrs,qgrs,xgrs,
#else
     &                  plamgr,pphigr,ugrs,vgrs,pgr,
     x                  tgrs,qgrs,xgrs,
#endif
#ifdef PSPLIT
     x                  ggt0,ggq0,ggu0,ggv0,
#else
     x                  gt0,gq0,gu0,gv0,
#endif
#ifdef DG
     x                  tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,
#endif
#ifdef DG3
     x                  gda,
#endif
#endif
#ifdef RAS
     x                  ras,lmx,cp,alhl,grav,rgas,
     x                  sig, sgb, prh, prj, hpk, fpk, ods, prns,
     x                  rasal, lm, krmin, krmax, nstrp,
     x                  ncrnd, rannum, afac, ufac,
#endif
#ifdef RASV2
     &                  ras,rgas, cp, grav, alhl,
     &                  sig,prj,sgb,rasal, rannum, dsfc,
     &                  pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,
#endif
#ifdef MUL_CLDTOP
     x                  xkt2,NCLDTOP,
#endif
#ifdef CLM_CWF
     x                  qugs, qvgs, cgs2, lwconv,
#endif
#ifdef SAS_DIAG
     +                  dcu,dcv,dct,dcq,dch,fcu,fcd,
     +                  deltb,delqb,delhb,cbmf,dlt,dlq,dlh,
#endif
#ifdef CONDENS
     +                  qcicnv,qrscnv,qcilrg,qrslrg,
#endif
#ifdef VIC
     x                  lat,count,idate)
#else
     x                  lat,count)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c
#ifdef CRAY_THREAD
cfpp$ expand(znlacm,accdia)
#endif
#ifndef RSM
#include <comfgrid.h>
#include <comfphys.h>
#include <comfver.h>
#ifdef ISOTOPE
#include <comfiso.h>
#endif
#include <comznl.h>
#include <comgpd.h>
       dimension idate(4)
#ifdef MP
#include <commpi.h>
#define LONF2S lonf2p_
#define LONF22S lonf22p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LONF22S lonf22_
#define LATG2S latg2_
#endif
#else
#include <rscomf.rerun.h>
#include <rscomltb.h>
#include <rscommap.h>
#include <rscomgpd.h>
#include <rscombgt.h>
#ifdef MP
#include <comrmpi.h>
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
#ifdef NFDR
#include <comnfdr.h>
      parameter(twopm=14.)
#endif
      parameter(cnwatt=-cal_*1.e4/60.,fv=rv_/rd_-1.,eps=rd_/rv_)
      parameter(epsm1=eps-1.,qmin=qmin_)
      parameter(iun=84)
#ifndef CRITSNOW
      parameter(critsnow=273.16)
#else
      parameter(critsnow=CRITSNOW)
#endif
#ifdef KEN
      dimension iptlat(nptken)
      dimension jptlat(nptken)
#endif
#ifndef RSM
c
#ifndef SCM
      dimension   pllamgr(LONF22S)
      dimension   plphigr(LONF22S)
#endif
      dimension      plgr(LONF22S)
      dimension     pslap(LONF22S)
      dimension     uugrs(LONF22S,levs_)
      dimension     vvgrs(LONF22S,levs_)
      dimension     tvgrs(LONF22S,levs_)
      dimension     rqgrs(LONF22S,levh_)
      dimension     xxgrs(LONF22S,levs_)
      dimension     ggt0(LONF22S,levs_)
      dimension     ggq0(LONF22S,levh_)
      dimension     ggu0(LONF22S,levs_)
      dimension     ggv0(LONF22S,levs_)
#else
#ifdef PSPLIT
      dimension     ggt0(LONF2S,levs_)
      dimension     ggq0(LONF2S,levh_)
      dimension     ggu0(LONF2S,levs_)
      dimension     ggv0(LONF2S,levs_)
#endif
#endif
c
c warning: gt0,gq0,gu0,gv0 may overlay tgrs,qgrs,ugrs,vgrs respectively.
c use local arrays to save the inputs
c
#ifndef SCM
      dimension   plamgr(LONF2S)
      dimension   pphigr(LONF2S)
#endif
      dimension     ugrs(LONF2S,levs_)
      dimension     vgrs(LONF2S,levs_)
      dimension      pgr(LONF2S)
      dimension     tgrs(LONF2S,levs_)
      dimension     qgrs(LONF2S,levh_)
      dimension     xgrs(LONF2S,levs_)
      dimension      gt0(LONF2S,levs_)
      dimension      gq0(LONF2S,levh_)
      dimension      gu0(LONF2S,levs_)
      dimension      gv0(LONF2S,levs_)
c
      dimension               slk(levs_)
c
      dimension kbot(LONF2S),
     *          ktop(LONF2S),
     *          kuo (LONF2S),
     *          dudt(LONF2S,levs_),
     *          dvdt(LONF2S,levs_),
     *          dtdt(LONF2S,levs_),
     *          dqdt(LONF2S,levh_),
     *          hsw(LONF2S,levs_),
     *          gflx(LONF2S),
     *          rain(LONF2S),
     *          rainc(LONF2S),
     *          rainl(LONF2S),
     *          rain1(LONF2S),
     *          rain2(LONF2S),
     *          evapc(LONF2S),
     *          work1(LONF2S),
     *          wind(LONF2S),
     *          w2(LONF2S,levs_),
     *          vvel(LONF2S,levs_)
c
#ifdef ISOTOPE
      dimension raini(LONF2S,2),
     &          rainci(LONF2S,2),
     &          rainli(LONF2S,2),
     &          rain1i(LONF2S,2)
      dimension condi(LONF2S,2),
     &          condci(LONF2S,2),
     &          condli(LONF2S,2),
     &          cond1i(LONF2S,2)
      dimension cond(LONF2S),
     *          condc(LONF2S),
     *          condl(LONF2S),
     *          cond1(LONF2S)
      dimension eqfo(LONF2S),eqfd(LONF2S)
      dimension revap(LONF2S,2)
      dimension qssi(LONF2S,2)
      dimension dqsfc1i1(LONF2S),dqsfc1i2(LONF2S)
      dimension tmpr(2),tmpc(2)
      dimension smsoil0(LONF2S,lsoil_)
      dimension rrsvi1(LONF2S,lsoil_)
      dimension rrsvi2(LONF2S,lsoil_)
#endif
c
#ifdef CO2PROG
** co2 **
      dimension dqsfc1i1(LONF2S)
#endif
      dimension  tstar(LONF2S), psexp(LONF2S),
     1          snowmt(LONF2S),factor(LONF2S),
     2          snowev(LONF2S),snowfl(LONF2S),
     3              fm(LONF2S),    fh(LONF2S),
     4              cd(LONF2S),   cdq(LONF2S),
     5             qss(LONF2S), radsl(LONF2S)
      dimension dusfcg(LONF2S),dvsfcg(LONF2S)
      dimension dusfc1(LONF2S),dvsfc1(LONF2S),
     *          dtsfc1(LONF2S),dqsfc1(LONF2S),
     *          dlwsf1(LONF2S),ulwsf1(LONF2S)
chybrid
      dimension sihyb(LONF2S,levp1_),cihyb(LONF2S,levp1_)
      dimension slhyb(LONF2S,levs_),clhyb(LONF2S,levs_)
      dimension delhyb(LONF2S,levs_),slkhyb(LONF2S,levs_)
      dimension prsi(levp1_),prsl(levs_)
#ifdef VIC
      dimension dswsf1(LONF2S)
#endif
#ifndef VIC
      dimension smsoil(LONF2S,lsoil_),
     *          stsoil(LONF2S,lsoil_)
#else
      dimension smsoil(LONF2S,lsoil_),
     *          stsoil(LONF2S,nsoil_)
#endif
      integer soiltyp,vegtype
      dimension soiltyp(LONF2S),
     *          sigmaf(LONF2S),
     *          vegtype(LONF2S)
      dimension rb(LONF2S),
     *          rhscnpy(LONF2S)
      dimension ai(LONF2S,lsoil_),
     *          bi(LONF2S,lsoil_),
     1          drain(LONF2S)
      dimension cci(LONF2S,lsoil_),
     *          rhsmc(LONF2S,lsoil_),
     *          runof(LONF2S)
      dimension zsoil(LONF2S,lsoil_),
     *          cld1d(LONF2S)
      dimension evap(LONF2S),
     *          hflx(LONF2S),
     *          rnet(LONF2S)
      dimension t850(LONF2S),
     *          pk(LONF2S),
     *          pkp(LONF2S),
     *          ep1d(LONF2S)
      dimension gamt(LONF2S),
     *          gamq(LONF2S),
     *          kpbl(LONF2S)
      dimension var(LONF2S),
     *          oc(LONF2S),
     1          oa4(LONF2S,4),
     *          ol4(LONF2S,4)
      dimension qcicnv(LONF2S,levs_),
     *           qrscnv(LONF2S,levs_)
      dimension ustar(LONF2S)
#ifdef HYDRO
      dimension hydrow(LONF2S,levs_)
      dimension evcn(LONF2S)
#else
#ifdef NOAHYDRO
      dimension evcn(LONF2S)
#endif
#endif
      dimension qu00(LONF2S),qv00(LONF2S)
#ifdef VIC
      dimension evp(LONF2S),heat(LONF2S),
     *          rho(LONF2S)
      real      ls, le
#endif
#ifdef NOALSM1
      dimension slptyp(LONF2S),
     1          slsoil(LONF2S,lsoil_)
      integer   slptyp
#endif
#ifdef MRGLSM
      dimension fm10   (LONF2S),fh2      (LONF2S)
      dimension oldprcp(LONF2S),oldsrflag(LONF2S)
#endif
#ifdef VICLSM1
      dimension expld(LONF2S,lsoil_),
     1          dphld(LONF2S,lsoil_),bubld(LONF2S,lsoil_),
     2          qrtld(LONF2S,lsoil_),bkdld(LONF2S,lsoil_),
     3          sldld(LONF2S,lsoil_),wcrld(LONF2S,lsoil_),
     4          wpwld(LONF2S,lsoil_),smrld(LONF2S,lsoil_),
     5          smxld(LONF2S,lsoil_),sicld(LONF2S,lsoil_),
     6          dpnld(LONF2S,nsoil_),sxnld(LONF2S,nsoil_),
     7          epnld(LONF2S,nsoil_),bbnld(LONF2S,nsoil_),
     8          apnld(LONF2S,nsoil_),btnld(LONF2S,nsoil_),
     9          gmnld(LONF2S,nsoil_)
      dimension vrtld(LONF2S,lsoil_),lstsnld(LONF2S)
      integer lstsnld
      real    kstld(LONF2S,lsoil_)
#endif
#ifdef VICLSM2
      dimension expld(LONF2S,lsoil_),
     1          dphld(LONF2S,lsoil_),bubld(LONF2S,lsoil_),
     2          qrtld(LONF2S,lsoil_),bkdld(LONF2S,lsoil_),
     3          sldld(LONF2S,lsoil_),wcrld(LONF2S,lsoil_),
     4          wpwld(LONF2S,lsoil_),smrld(LONF2S,lsoil_),
     5          smxld(LONF2S,lsoil_)
     6          dpnld(LONF2S,nsoil_),sxnld(LONF2S,nsoil_),
     7          epnld(LONF2S,nsoil_),bbnld(LONF2S,nsoil_),
     8          apnld(LONF2S,nsoil_),btnld(LONF2S,nsoil_),
     9          gmnld(LONF2S,nsoil_)
      real    kstld(LONF2S,lsoil_)

c
      dimension flaild(LONF2S,msub_)
      dimension vfrld(LONF2S,msub_), vtyld(LONF2S,msub_)
      dimension cnpld(LONF2S,msub_), snold(LONF2S,msub_)
      dimension csnld(LONF2S,msub_), rsnld(LONF2S,msub_)
      dimension tsfld(LONF2S,msub_), tpkld(LONF2S,msub_)
      dimension sfwld(LONF2S,msub_), pkwld(LONF2S,msub_)
c
      dimension sicld(LONF2S,lsoil_,msub_)
      dimension smcld(LONF2S,lsoil_,msub_)
      dimension stcld(LONF2S,nsoil_,msub_)
c
      dimension vrtld(LONF2S,lsoil_,msub_)
      dimension lstsnld(LONF2S,msub_), nvegld(LONF2S)
      integer lstsnld, nvegld
#endif
#ifdef CONDENS
      dimension qcilrg(LONF2S,levs_),qrslrg(LONF2S,levs_)
      dimension cldwp(LONF2S,levs_),cldip(LONF2S,levs_)
#endif

#ifdef ZM2
      common /zm2/
     1 tpzm(LONF2S,levs_,LATG2S)
     2,tp1zm(LONF2S,levs_,LATG2S)
     3,qpzm(LONF2S,levs_,LATG2S)
     4,qp1zm(LONF2S,levs_,LATG2S)
     5,idpinizm(LATG2S)
      data idpinizm /LATG2S*0/
#endif

#ifdef CLD2
      common /cloud2/
     1 tp(LONF2S,levs_,LATG2S)
     2,tp1(LONF2S,levs_,LATG2S) 
     3,qp(LONF2S,levs_,LATG2S)
     4,qp1(LONF2S,levs_,LATG2S)
     5,psp(LONF2S,LATG2S)
     6,psp1(LONF2S,LATG2S)
     7,idpini(LATG2S)
       data idpini /LATG2S*0/
#endif
      dimension tovzer(levs_)
#ifdef DG3
#include <comgda.h>
      dimension gda(nwgda,kdgda)
#endif
#ifdef MUL_CLDTOP
      dimension xkt2(LONF2S*NCLDTOP)
#endif
#ifdef SAS_DIAG
      dimension dcu(LONF2S,levs_),dcv(LONF2S,levs_),dct(LONF2S,levs_),
     +          dcq(LONF2S,levs_),dch(LONF2S,levs_),
     +          fcu(LONF2S,levs_),fcd(LONF2S,levs_)
      dimension deltb(LONF2S),delqb(LONF2S),delhb(LONF2S),cbmf(LONF2S)
      dimension dlt(LONF2S,levs_),dlq(LONF2S,levs_),dlh(LONF2S,levs_)
#endif
#ifdef RAS
      logical ras
      dimension sig(lmx+1), prj(lmx+1), prh(lmx),   fpk(lmx), hpk(lmx)
     *,         sgb(lmx),   ods(lmx),   rasal(lmx), prns(lmx/2)
     *,         rannum(200)
      real    clw(LONF2S,levs_)  ! detrained cloud water at cloud top [kg/kg]
      real    cmf(LONF2S,levs_)  ! cumulus mass flux []
      real    detr(LONF2S,levs_) ! mass detrainment [-]

#ifdef ISOTOPE
      dimension gt00(levs_),gq00(levh_),gu00(levs_),gv00(levs_)
#else
      dimension gt00(levs_),gq00(levs_),gu00(levs_),gv00(levs_)
#endif
      dimension clw0(levs_),cmf0(levs_),detr0(levs_)
      dimension si0(levs_+1),sl0(levs_),del0(levs_)
#endif
#ifdef RASV2
      logical ras
      real sig(levs_+1), prj(levs_+1), sgb(levs_), rasal(levs_)
     *,    rannum(ncrnd), dsfc(LONF2S),pdd
      integer krmin, krmax, kfmax, ncrnd, mct, kctop(mct+1)
      dimension gt00(levs_),gq00(levs_),gu00(levs_),gv00(levs_)
      dimension clw0(levs_),clt0(levs_)
      dimension clw(LONF2S,levs_),clt(LONF2S,levs_)
#endif
#ifdef KF2
      dimension clw(LONF2S,levs_),cli(LONF2S,levs_)
#endif
#ifdef ZM2
      real     clw(LONF2S,levs_)  ! detrained cloud water [kg/kg]
      real     detr(LONF2S,levs_) ! mass detrainment      [-]
#endif

      logical ldiag
c
#ifndef NOAHYDRO
      parameter (hfus=hfus_,rhoh2o=1000.,epsq=2.e-12)
#else
      parameter (hfus=hfus_,hvap=hvap_,rhoh2o=1000.,epsq=2.e-12)
#endif
c
#ifdef CLDT93
      real    hdf(LONF2S,levs_)
      real    qdf(LONF2S,levs_)
#endif
#ifdef DBG
      real      qck1(LONF2S,ncldg_+1)
      real      qck2(LONF2S,ncldg_+1)
      real      rck1(LONF2S),rck2(LONF2S)
      real      qsum1, qsum2, qsum3
      parameter (ucnv1=1000.0/g_ ,ucnv2=1000.0)
#endif
#ifdef NFDR
      real rh2m(LONF2S)
#endif
c
#ifdef CO2PROG
** co2 **
      common /co2/co2fl(LONF2S,LATG2S)
#endif
#ifdef RIVER
#include <comfriv.h>
#endif
c
      ncloud = ncldg_
      ntrace = ntrac_
      icloud = ntrac_ + 1
      kcloud = ntrac_*levs_+1
      ktrace = levs_+1
      ktotal = levs_*ntotal_
c
      if(ncloud.le.0) then
         nvdiff = ntrac_
      elseif(ncloud.le.2) then
         nvdiff = ntrac_ + 1
      elseif(ncloud.le.5) then
         nvdiff = ntrac_ + 2
      endif
#ifndef RSM
#define LONS2 lons2
#else
#ifdef MP
#ifdef RSMVECTORIZE
#define LONS2 igrd12p_*jgrd12p_
#else
#define LONS2 lonlen(mype)*2
#endif
#else
#define LONS2 igrd12_
#endif
#endif
#ifndef RSM
c
c   reserve the input without change
c
#ifdef CLM_CWF
      do k=1,levs_
        do j=1,lons2
          qugs(j,lat,k)=uugrs(j,k)*rqgrs(j,k)
          qvgs(j,lat,k)=vvgrs(j,k)*rqgrs(j,k)
        enddo
      enddo
#endif


      do j=1,LONS2
#ifndef SCM
        plamgr(j)=pllamgr(j)
        pphigr(j)=plphigr(j)
#endif
        pgr   (j)=plgr   (j)
      enddo
      do k=1,levs_
        do j=1,LONS2
          ugrs(j,k)=uugrs(j,k)
          vgrs(j,k)=vvgrs(j,k)
          tgrs(j,k)=tvgrs(j,k)
          xgrs(j,k)=xxgrs(j,k)
#ifdef ISOTOPE
        enddo
      enddo
      do k=1,levh_
        do j=1,LONS2
#endif
          qgrs(j,k)=rqgrs(j,k)
        enddo
      enddo
      if(ntotal_.gt.1) then
        if (kcloud.gt.ktrace) then
          do k=ktrace,kcloud-1
            do j=1,LONS2
              qgrs(j,k)=rqgrs(j,k)  ! tracer
            enddo
          enddo
        end if
        if (levh_.gt.kcloud) then
          do k=kcloud,levh_
            do j=1,LONS2
#ifndef NISLQ
#ifdef CLDT93
              qgrs(j,k)=max(rqgrs(j,k),0.0)  ! cloud water ...
#else
              qgrs(j,k)=max(rqgrs(j,k),qmin) ! cloud water ...
#endif
#else  /* NISLQ */
#ifdef CLDT93
              qgrs(j,k)=max(rqgrs(j,k),0.D0)  ! cloud water ...
#else
              qgrs(j,k)=max(rqgrs(j,k),qmin) ! cloud water ...
#endif
#endif
            enddo
          enddo
        end if
      endif

#ifndef O3C
      if(ntrace.ge.2.and.thour.eq.0.) then
#ifdef DBG
#ifndef MP
       if(lat.eq.LATG2S)
     1 call rmaxmin(qgrs(1,ktrace),
     1              LONS2,1,1,1,'trace before')
#endif
#endif
        do k=1,levs_
          do j=1,LONS2
            qgrs(j,k+ktrace-1)=ozon(j,k,lat)
          enddo
        enddo
#ifdef DBG
#ifndef MP
       if(lat.eq.LATG2S)
     1 call rmaxmin(qgrs(1,ktrace),
     1              LONS2,1,1,1,'trace after')
#endif
#endif
      endif

#endif
#endif
c
#ifdef ISOTOPE
*** isotope correction 1 ***
      do k=1,levs_
        do j=1,LONS2
          qgrs(j,k+levs_)=0.D0 !! dummy
          if (qgrs(j,k).le.0.D0) then
            qgrs(j,k)=0.D0
            qgrs(j,k+2*levs_)=0.D0
            qgrs(j,k+3*levs_)=0.D0
          else
            qgrs(j,k+2*levs_)=max(qgrs(j,k+2*levs_),0.D0)  ! tracer
            qgrs(j,k+3*levs_)=max(qgrs(j,k+3*levs_),0.D0)  ! tracer
          endif
c
#ifndef NOFRAC
          if ((qgrs(j,k).le.1.D-5).and.(qgrs(j,k).gt.0.D0).and.
     $           ((qgrs(j,k+2*levs_)/qgrs(j,k).gt.1.05D0).or.
     $            (qgrs(j,k+3*levs_)/qgrs(j,k).gt.1.40D0))) then
            alpha=1.-1./(-log10(qgrs(j,k)/1.D-5)+1)**2.
            qgrs(j,k+2*levs_)=qgrs(j,k+2*levs_)+
     $           (qgrs(j,k)-qgrs(j,k+2*levs_))*alpha
            qgrs(j,k+3*levs_)=qgrs(j,k+3*levs_)+
     $           (qgrs(j,k)-qgrs(j,k+3*levs_))*alpha
c            print'(a6,4i4,3e15.7)','gbCHK+',mype,j,k,lat,
c     $            qgrs(j,k),
c     $            qgrs(j,k+2*levs_)/qgrs(j,k),
c     $            qgrs(j,k+3*levs_)/qgrs(j,k)
          endif
#endif
c          if ((mype.eq.58).and.(j.eq.7).and.(lat.eq.1)) then
c             print'(a6,4i3,5e15.7)','BEFPHS',mype,j,k,lat,
c     $            qgrs(j,k),qgrs(j,k+2*levs_),
c     $            qgrs(j,k+2*levs_)/qgrs(j,k),
c     $            qgrs(j,k+3*levs_),qgrs(j,k+3*levs_)/qgrs(j,k)
c          endif
        enddo
      enddo
#endif
c
      sigshc=0.7
      levshc=0
      p850 = 850. * .1
#ifdef RSM
      do j=1,LATG2S
        rbs2 = 1.0
      enddo
      jcapr= 126.*106000./ rdelx
      dxmeter = rdelx
#else
#ifdef SCM
      jcapset = 126
      dxmeter = 12600000./jcapset
#else
      dxmeter = 12600000./jcap_
#endif
#endif
c
c  frain is the factor for centered difference scheme correction of rain
c  amount.
c
      frain = .5
      if(inistp.eq.1) frain = 1.
      do 205 k=1,levs_
        slk(k)=sl(k)**(rd_/cp_)
205   tovzer(k)=0. e 0
      msta=0
      ldiag=.false.
c
c  transfer vegetation fraction from global to local location
c  as well as vegetation type and soil type
c
      do j = 1, LONS2
#ifndef R2
#ifdef OSULSM2
        soiltyp(j)= stype(j,lat)+.5
#ifdef NOALSM1
        sigmaf(j) = max(vfrac(j,lat),0.)
#else
#ifdef NCAR_EDIR
        sigmaf(j) = max(vfrac(j,lat),0.)
#else
        sigmaf(j) = max(vfrac(j,lat),.3)
#endif
#endif  /* ifndef NOALSM1 */
#endif  /* ifdef OSULSM2 */
#ifdef NOA
        soiltyp(j)= stype(j,lat)+.5
#ifdef NOALSM1
        sigmaf(j) = max(vfrac(j,lat),0.)
#else
#ifdef NCAR_EDIR
        sigmaf(j) = max(vfrac(j,lat),0.)
#else
        sigmaf(j) = max(vfrac(j,lat),.3)
#endif
#endif  /* ifndef NOALSM1 */
#endif  /* ifdef NOA */
#ifdef VIC
        sigmaf(j) = max(vfrac(j,lat),0.)
#endif
#ifdef OSULSM1
        soiltyp(j) = 7
        sigmaf(j) = .7
#else
        vegtype(j)= vtype(j,lat)+.5
        fm(j) = ffmm(j,lat)
        fh(j) = ffhh(j,lat)
        ustar(j) = uustar(j,lat)
#endif
#ifdef OSULSM2
        if(slmsk(j,lat).eq.2.) then
#ifdef STATSGO_SOIL
          vegtype(j) = 12
          soiltyp(j) = 16
#else
          vegtype(j) = 13
          soiltyp(j) = 9
#endif
        endif
#endif
#endif   /* ifndef R2  */
#ifdef R2
        soiltyp(j) = 7
        sigmaf(j) = .7
#ifndef OSULSM1
        vegtype(j)= vtype(j,lat)+.5
        fm(j) = ffmm(j,lat)
        fh(j) = ffhh(j,lat)
        ustar(j) = uustar(j,lat)
        if(slmsk(j,lat).eq.2.) then
          vegtype(j) = 13
          soiltyp(j) = 9
        endif
#endif
#endif
#ifdef NOALSM1
        slptyp(j) = slope(j,lat)+.5
        if(slmsk(j,lat).eq.2.) then
          slptyp(j) = 0
        endif
#endif
      enddo
c
      sl1=sl(1)
      factx=1000.e0*(deltim*del(1)/g_)
      cstrn2=1. e 0/(.75 e 0*150. e 0)
      dt2=deltim*2. e 0
      dtf=frain*dt2
#ifndef RSM
      dtf=dtf*count
#endif
c
c  transfer soil moisture and temperature from global to local variables
c
      do k = 1, lsoil_
        do j = 1, LONS2
          smsoil(j,k) = smc(j,lat,k)
#ifndef VIC
          stsoil(j,k) = stc(j,lat,k)
#endif
#ifdef NOALSM1
          slsoil(j,k) = slc(j,lat,k)
#endif
#ifdef VIC
          expld(j,k) = expt(j,lat,k)
          kstld(j,k) = kest(j,lat,k)
          dphld(j,k) = dph(j,lat,k)
          bubld(j,k) = bub(j,lat,k)
          qrtld(j,k) = qrt(j,lat,k)
          bkdld(j,k) = bkd(j,lat,k)
          sldld(j,k) = sld(j,lat,k)
          wcrld(j,k) = wcr(j,lat,k)
          wpwld(j,k) = wpw(j,lat,k)
          smrld(j,k) = smr(j,lat,k)
          smxld(j,k) = smx(j,lat,k)
          vrtld(j,k) = vroot(j,lat,k)
#endif
#ifdef VICLSM1
          sicld(j,k) = sic(j,lat,k)
#endif
#ifdef ISOTOPE
          smsoil0(j,k)=smc(j,lat,k)
          rrsvi1(j,k)=rresv1(j,lat,k)
          rrsvi2(j,k)=rresv2(j,lat,k)
#endif
        enddo
      enddo
c
#ifdef VICLSM2
      do j = 1, LONS2
         nvegld(j) = nveg(j,lat)
         do nv = 1, nvegld(j)
            flaild(j,nv) = flai(j,lat,nv)
            vfrld(j,nv)  = mvfr(j,lat,nv)      ! tile coverage
            vtyld(j,nv)  = mvty(j,lat,nv)
            cnpld(j,nv)  = mcnp(j,lat,nv)
            snold(j,nv)  = msno(j,lat,nv)      ! tile snow depth
            csnld(j,nv)  = csno(j,lat,nv)      ! canopy snow
            rsnld(j,nv)  = rsno(j,lat,nv)
            tsfld(j,nv)  = tsf(j,lat,nv)
            tpkld(j,nv)  = tpk(j,lat,nv)
            sfwld(j,nv)  = sfw(j,lat,nv)
            pkwld(j,nv)  = pkw(j,lat,nv)
            lstsnld(j,nv)= lstsn(j,lat,nv)
         end do
      enddo
      do k = 1, lsoil_
         do j = 1, LONS2
            do nv = 1, nvegld(j)
               vrtld(j,k,nv)=vroot(j,lat,(nv-1)*lsoil_+k)
               smcld(j,k,nv)=msmc(j,lat,(nv-1)*lsoil_+k)
               sicld(j,k,nv)=msic(j,lat,(nv-1)*lsoil_+k)
            enddo
         enddo
      end do
#endif
c
#ifdef VIC
      do k = 1, nsoil_
         do j = 1, LONS2
            stsoil(j,k)= stc(j,lat,k)
            dpnld(j,k) = dphn(j,lat,k)
            sxnld(j,k) = smxn(j,lat,k)
            epnld(j,k) = expn(j,lat,k)
            bbnld(j,k) = bubn(j,lat,k)
            apnld(j,k) = alpn(j,lat,k)
            btnld(j,k) = betn(j,lat,k)
            gmnld(j,k) = gamn(j,lat,k)
         enddo
      enddo
#endif
#ifdef VICLSM2
      do k = 1, nsoil_
         do j = 1, LONS2
            do nv = 1, nvegld(j)
               stcld(j,k,nv)=mstc(j,lat,(nv-1)*nsoil_+k)
            end do
         end do
      end do
#endif
#ifdef VICLSM1
      do j = 1, LONS2
        lstsnld(j) = -1
      enddo
      do j = 1, LONS2
        if(slmsk(j,lat).eq.1) then
           if(lstsn(j,lat).gt.1.e10) then
             lstsnld(j) = -1
           else
             lstsnld(j) = lstsn(j,lat)
           endif
        endif
      end do
#endif
c
      do k=1,levs_
        do j=1,LONS2
          dudt(j,k)=0. e 0
          dvdt(j,k)=0. e 0
          dtdt(j,k)=0. e 0
        enddo
      enddo
      do k=1,levh_
        do j=1,LONS2
          dqdt(j,k)=0. e 0
        enddo
      enddo
#ifndef RSM
      rcl=rbs2(lat)
#else
      rcl=1.
#endif
      rcs = sqrt(rcl)
c
c  get dry tep.
c
#ifdef DBG
#ifndef MP
      if(lat.eq.LATG2S) 
     1    call rmaxmin(tgrs(1,1),LONS2,
     1    1,1,1,'tgrs virtual')
#endif
#endif
      do 260 k=1,levs_
      do 260 j=1,LONS2
        w2(j,k)=1.+fv*qgrs(j,k)
        work1(j) = 0.
        if(ncloud.gt.0) then
           do ic = icloud,ntotal_
              kc = (ic-1)*levs_ + k
#ifndef NISLQ
              work1(j) = max(qgrs(j,kc),qmin) + work1(j)
#else
              work1(j) = max(qgrs(j,kc),qmin) + work1(j)
#endif
           enddo
           w2(j,k)=w2(j,k) - work1(j)
        endif
  260 continue
      do 270 k=1,levs_
      do 270 j=1,LONS2
        tgrs(j,k)=tgrs(j,k)/w2(j,k)
  270 continue
#ifdef DBG
#ifndef MP
      if(lat.eq.latg2_) 
     1    call rmaxmin(tgrs(1,1),LONS2,
     1                 1,1,1,'tgrs dry')
#endif
#endif
#ifdef DG
      do 271 k=1,levs_
        j=ismax(LONS2,tgrs(1,k),1)
        if(tgrs(j,k).gt.tgmxl) then
          tgmxl=tgrs(j,k)
          igmxl=j
          kgmxl=k
        endif
        j=ismin(LONS2,tgrs(1,k),1)
        if(tgrs(j,k).lt.tgmnl) then
          tgmnl=tgrs(j,k)
          igmnl=j
          kgmnl=k
        endif
#endif
  271 continue
c
c  get surface pressure.
c
      do 280 j=1,LONS2
      psexp(j)= exp(pgr(j))
      factor(j)=factx*psexp(j)
      psurf(j,lat)=psexp(j)
      psmean(j,lat)=psmean(j,lat)+psexp(j)*dtf
  280 continue
#ifdef HYBRID
      do j=1,LONS2
        do k=1,levp1_
          sihyb(j,k)=ak5(levp1_-k+1)/psexp(j)+bk5(levp1_-k+1)
          cihyb(j,k)=1.-sihyb(j,k)
        enddo
        do k=1,levs_
          slhyb(j,k)=0.5*(sihyb(j,k)+sihyb(j,k+1))
          delhyb(j,k)=sihyb(j,k)-sihyb(j,k+1)
          clhyb(j,k)=1.-slhyb(j,k)
        enddo
      enddo
#else
      do j=1,LONS2
        do k=1,levp1_
          sihyb(j,k)=si(k)
          cihyb(j,k)=ci(k)
        enddo
        do k=1,levs_
          slhyb(j,k)=sl(k)
          delhyb(j,k)=del(k)
          clhyb(j,k)=cl(k)
        enddo
      enddo
#endif
      do j=1,LONS2
        do k=1,levs_
          slkhyb(j,k)=slhyb(j,k)**(rd_/cp_)
        enddo
      enddo
      do k=1,levs_
        do j=1,LONS2
        if((slhyb(j,k).gt.sigshc).and.(k.gt.levshc-1)) levshc=k+1
        enddo
      enddo
c
c  initialize dtdt with heating rate from dcyc2 and get radsl for progtm
c
#ifdef LFC
      if(1.ne.1) then
#endif
        call dcyc2(LONS2,LONF2S,
     1             levs_,
#ifndef RSM
     &             solhr,slag,sinlab(1,lat),coslab(1,lat),
     &             sdec,cdec,
     &             xlon(1,lat),coszen(1,lat),
#else
     &             solhr,slag,sinlar(1,lat),coslar(1,lat),
     &             sdec,cdec,
     &             rlon(1,lat),coszer(1,lat),
#endif
     &             sfcdlw(1,lat),sfcnsw(1,lat),
#ifdef VIC
     &             sfcdsw(1,lat),dswsf1,
#endif
     &             tgrs(1,1),
     &             tsea(1,lat),tsflw(1,lat),swh(1,1,lat),hlw(1,1,lat),
     &             dlwsf1,ulwsf1,radsl,dtdt)
#ifdef DBG
#ifndef MP
      if(lat.eq.latg2_) 
     1       call rmaxmin(tgrs(1,1),LONS2,
     1                    1,1,1,'dcyc2')
#endif
#endif
#ifdef LFC
      else
        do i=1,LONS2
          radsl(i)=sfcnsw(i,lat)+sfcdlw(i,lat)
          ulwsf1(i)=sbc_*tsea(i,lat)**4
        enddo
        do k=1,levs_
          do i=1,LONS2
            dtdt(i,k)=dtdt(i,k)+swh(i,k,lat)+hlw(i,k,lat)
          enddo
        enddo
      endif
#endif
c
      do 285 j=1,LONS2
      dlwsfc(j,lat)=dlwsfc(j,lat)+dlwsf1(j)*dtf
      ulwsfc(j,lat)=ulwsfc(j,lat)+ulwsf1(j)*dtf
285   continue
c
      do 286 k=1,levs_
      do 286 j=1,LONS2
      hsw(j,k)=dtdt(j,k)-hlw(j,k,lat)
286   continue
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            hsw,psexp, dtf,nmthsw,lat)
      call znlacm(LONS2,LONF2S,
     1            hlw(1,1,lat),psexp, 
     1            dtf,nmthlw,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            hsw, dtf,kdthsw,gda)
      call accdia(LONS2,LONF2S,
     1            hlw(1,1,lat), dtf,kdthlw,gda)
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'----- gbphys -----  lat=',lat
#ifdef MP
      endif
#endif
      level=1
#ifdef MP
      if(mype.eq.master) then
#endif
c      numlevs=levs_
      numlevs=1
      call maxmin(ugrs(1,level),LONS2,
     1            levs_,1,numlevs,'ugrs')
      call maxmin(vgrs(1,level),LONS2,
     1            levs_,1,numlevs,'vgrs')
      call maxmin(tgrs(1,level),LONS2,
     1            levs_,1,numlevs,'tgrs')
      call maxmin(qgrs(1,level),LONS2,
     1            levs_,1,numlevs,'qgrs')
#ifdef ISOTOPE
      call maxmin(qgrs(1,level+2*levs_),LONS2,
     1            levs_,1,numlevs,'qgrs2')
      call maxmin(qgrs(1,level+3*levs_),LONS2,
     1            levs_,1,numlevs,'qgrs3')
#endif
      call maxmin(hsw(1,level),LONS2,
     1            levs_,1,numlevs,'hsw')
      call maxmin(hlw(1,level,lat),LONS2,
     1            levs_,1,numlevs,'hlw')
      call maxmin(sheleg(1,lat),LONS2,
     1            1,1,1,'sheleg')
      call maxmin(tsea(1,lat),LONS2,
     1            1,1,1,'tsea')
      call maxmin(smsoil(1,level),LONS2,
     1            1,1,1,'smsoil')
      call maxmin(stsoil(1,level),LONS2,
     1            1,1,1,'stsoil')
      call maxmin(canopy(1,lat),LONS2,
     1            1,1,1,'canopy')
      call maxmin(zorl(1,lat),LONS2,
     1            1,1,1,'zorl')
#ifdef OSULSM1
      call maxmin(plantr(1,lat),LONS2,
     1            1,1,1,'plantr')
#endif
      call maxmin(tg3(1,lat),LONS2,
     1            1,LONS2,1,'tg3')
      call maxmin(slmsk(1,lat),LONS2,
     1            1,LONS2,1,'slmsk')
#ifdef MP
      endif
#endif
#endif

#ifdef MRGLSM
c
c Run merged LSM: sfcdif -> sfcdrv -> sfcdiag
c
c assign prcp and srflag from global to local array
c
      do j=1,LONS2
        oldprcp(j)   = prcp(j,lat)
        oldsrflag(j) = srflag(j,lat)
      enddo
c
#ifdef VIC
      call vday(idate(4),idate(2),idate(3),idate(1),thour,
     1          jyy,jmm,jdd,jhh,dummy)
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(tstar,LONS2,
     1            1,1,1,'tstar bef sfcdrv')
      call maxmin(qss,LONS2,
     1            1,1,1,'qss bef sfcdrv')
#ifdef MP
      endif
#endif
#ifdef MP
      if(mype.eq.master) then
#endif
      print*,' before sfcdrv lat = ',lat
      print*,LONS2,LONF2S,
     1 lsoil_,psexp(1),
     1 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     2 sheleg(1,lat),tsea(1,lat),qss(1),
     3 smsoil(1,1),stsoil(1,1),evapc(1),soiltyp(1),sigmaf(1),
     4 vegtype(1),canopy(1,lat),
     5 dlwsf1(1),radsl(1),
     6 deltim,zorl(1,lat),tg3(1,lat),
     7 sl1,slk(1),slmsk(1,lat),
     & wind(1),evcn(1)
#ifdef MP
      endif
#endif
#ifdef VIC
#ifdef MP
      if(mype.eq.master) then
#endif
       print *,'for viclsm variables'
       print *,nsoil_, jmm,
     + binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),
     + expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,smrld,
     + smxld,sicld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,
     + flai(1,lat),vrtld,lstsnld,
     + silz(1,lat),snwz(1,lat),csno(1,lat),rsno(1,lat),
     + tsf(1,lat),tpk(1,lat),sfw(1,lat),pkw(1,lat)
#ifdef MP
      endif
#endif
#endif
#endif
c
c compute exchange coeff (CD, CDQ, FFMM, FFHH)
c
      call sfcdif(LONS2,LONF2S,
     1            lsoil_,psexp,
     +            ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     +            tsea(1,lat), zorl(1,lat),
     +            cd,cdq,rb,
     +            rcl,sl1,slk(1),slmsk(1,lat),inistp,lat,
#ifndef RASV2
     +            fm,fh,ustar,wind,
#else
     +            fm,fh,ustar,wind,dsfc,
#endif
#ifdef VIC
     +            fm10,fh2,rho)
#else
     +            fm10,fh2)
#endif
c
c compute sfc energy/water balance
c
       call sfcdrv(LONS2,LONF2S,
     1             lsoil_,
     1             psexp,tgrs(1,1),qgrs(1,1),
     2             sheleg(1,lat),tsea(1,lat),qss,
     3             smsoil,stsoil,evapc,soiltyp,sigmaf,vegtype,
     4             canopy(1,lat),
     4             dlwsf1,radsl,snowmt,snowev,
     5             deltim,zorl(1,lat),tg3(1,lat),
     7             gflx,zsoil,cd,cdq,
     8             rcl,sl1,slk(1),slmsk(1,lat),inistp,lat,
#ifndef NOAHYDRO
     +             drain,evap,hflx,ep1d,wind,
#else
     +             drain,evap,hflx,ep1d,wind,evcn,
#endif
#ifdef NOALSM1
     +             snwdph(1,lat),slsoil,snoalb(1,lat),
     +             slptyp,shdmin(1,lat),shdmax(1,lat),
#endif
#ifdef VICLSM1
     + nsoil_, jmm,dswsf1,
     + binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),
     + expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,smrld,
     + smxld,sicld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,
     + flai(1,lat),vrtld,lstsnld,
     + silz(1,lat),snwz(1,lat),csno(1,lat),rsno(1,lat),
     + tsf(1,lat),tpk(1,lat),sfw(1,lat),pkw(1,lat),
#endif
#ifdef VICLSM2
     + nsoil_,msub_,jmm,dswsf1,
     + binf(1,lat),ds(1,lat),dsm(1,lat),ws(1,lat),cef(1,lat),
     + expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,wpwld,
     + smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,btnld,gmnld,
     + nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,rsnld,
     + tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,smcld,sicld,stcld,
     + silz(1,lat),snwz(1,lat),
#endif
     + snowfl, runof,
     + oldprcp,oldsrflag)
c
c update surafce layer properties
c
      call sfcdiag(LONS2,LONF2S,
     1             lsoil_,
     1 psexp,ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),tsea(1,lat),
#ifdef ISOTOPE
     $ qgrs(1,1+2*levs_),qgrs(1,1+3*levs_),
#endif
     4 radsl,rnet,
     6 f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),q2m(1,lat),
#ifdef ISOTOPE
     $ q2mi1(1,lat),q2mi2(1,lat),
#endif
     8 rcl,slk(1),slmsk(1,lat),inistp,lat,
     9 qss,evap,fm,fh,fm10,fh2)
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print*,' after sfcdrv lat = ',lat
      print*,LONS2,LONF2S,
     1 lsoil_,psexp(1),
     1 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     2 sheleg(1,lat),tsea(1,lat),qss(1),
     3 smsoil(1,1),stsoil(1,1),evapc(1),soiltyp(1),sigmaf(1),
     1 vegtype(1),canopy(1,lat),
     4 dlwsf1(1),radsl(1),snowmt(1),snowev(1),
     5 deltim,zorl(1,lat),tg3(1,lat),
     6 gflx(1),f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),
     & q2m(1,lat),
     7 (zsoil(1,k),k=1,lsoil_),cd(1),cdq(1),
     7 rcl,sl1,slk(1),slmsk(1,lat),lat,
     9 drain(1),evap(1),hflx(1),rnet(1),ep1d(1),fm(1),fh(1),ustar(1),
     & wind(1),evcn(1)
#ifdef MP
      endif
#endif
#endif
#else
c
c  Run 2-step LSM (progtm99 or progtm)
c
#ifndef OSULSM1
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print*,' before progtm lat = ',lat
      print*,LONS2,LONF2S,
     1 lsoil_,psexp(1),
     1 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     2 sheleg(1,lat),tsea(1,lat),qss(1),
     3 smsoil(1,1),stsoil(1,1),evapc(1),soiltyp(1),sigmaf(1),
     4 vegtype(1),canopy(1,lat),
     5 dlwsf1(1),radsl(1),
     6 deltim,zorl(1,lat),tg3(1,lat),
     7 sl1,slk(1),slmsk(1,lat)
#ifdef MP
      endif
#endif
#endif
      call progtm99(LONS2,LONF2S,
     1 lsoil_,psexp,
     1 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     2 sheleg(1,lat),tsea(1,lat),qss,
     3 smsoil,stsoil,evapc,soiltyp,sigmaf,vegtype,canopy(1,lat),
     4 dlwsf1,radsl,snowmt,snowev,
     5 deltim,zorl(1,lat),tg3(1,lat),
     6 gflx,f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),q2m(1,lat),
     7 zsoil,cd,cdq,rb,rhscnpy,rhsmc,ai,bi,cci,
#ifndef VIC
     8 rcl,sl1,slk(1),slmsk(1,lat),inistp,lat,
#else
     8 rcl,sl1,slk(1),slmsk(1,lat),inistp,rho,lat,
#endif
#ifdef HYDRO
#ifndef RASV2
     9 drain,evap,hflx,rnet,ep1d,fm,fh,ustar,wind,evcn,hydrow)
#else
     9 drain,evap,hflx,rnet,ep1d,fm,fh,ustar,wind,dsfc,evcn,hydrow)
#endif
#else
#ifndef RASV2
     9 drain,evap,hflx,rnet,ep1d,fm,fh,ustar,wind)
#else
     9 drain,evap,hflx,rnet,ep1d,fm,fh,ustar,wind,dsfc)
#endif
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print*,' after progtm99 lat = ',lat
      print*,LONS2,LONF2S,
     1 lsoil_,psexp(1),
     1 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),
     2 sheleg(1,lat),tsea(1,lat),qss(1),
     3 smsoil(1,1),stsoil(1,1),evapc(1),soiltyp(1),sigmaf(1),
     1 vegtype(1),canopy(1,lat),
     4 dlwsf1(1),radsl(1),snowmt(1),snowev(1),
     5 deltim,zorl(1,lat),tg3(1,lat),
     6 gflx(1),f10m(1,lat),u10m(1,lat),v10m(1,lat),t2m(1,lat),
     & q2m(1,lat),
     7 (zsoil(1,k),k=1,lsoil_),cd(1),cdq(1),
     7 rcl,sl1,slk(1),slmsk(1,lat),lat,
     9 drain(1),evap(1),hflx(1),rnet(1),ep1d(1),fm(1),fh(1),ustar(1),
     & wind(1)
#ifdef MP
      endif
#endif
#endif
#else
      call progtm(LONS2,LONF2S,
     *            lsoil_,psexp,
     *            ugrs(1,1),vgrs(1,1),
     *            tgrs(1,1),qgrs(1,1),
     *            sheleg(1,lat),tsea(1,lat),qss,
     *            smsoil,stsoil,evapc,soiltyp,sigmaf,canopy(1,lat),
     *            radsl,snowmt,snowev,
     *            deltim,zorl(1,lat),plantr(1,lat),tg3(1,lat),
     *            gflx,f10m(1,lat),u10m(1,lat),v10m(1,lat),
     *            t2m(1,lat),q2m(1,lat),
     *            zsoil,cd,cdq,rb,rhscnpy,rhsmc,ai,bi,cci,
#ifndef VIC
     *            rcl,sl1,slk(1),slmsk(1,lat),inistp,lat,
#else
     *            rcl,sl1,slk(1),slmsk(1,lat),inistp,lat,rho,
#endif
#ifndef RASV2
     *            drain,evap,hflx,rnet,ep1d,cowave,fm,fh,wind)
#else
     *            drain,evap,hflx,rnet,ep1d,cowave,fm,fh,wind,dsfc)
#endif
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(tsea(1,lat),LONS2,
     1            1,1,1,'tsea after progtm')
#ifdef MP
      endif
#endif
#endif
#endif
c
#ifdef RASV2
      do j=1,LONS2
        dsfc(j)=0.
      enddo
#endif
#ifdef MRGLSM
      do j=1,LONS2
        snowfall(j,lat) = snowfall(j,lat) + snowfl(j)
      enddo
#endif
      do j=1,LONS2
#ifdef HYDRO
        snowevap(j,lat) = snowevap(j,lat) + snowev(j)*dtf
#else
#ifdef NOAHYDRO
        if (slmsk(j,lat).eq.1) 
     &      snowevap(j,lat) = snowevap(j,lat) 
     &      + snowev(j)*deltim*(hvap+hfus)
        if (slmsk(j,lat).eq.2) 
     &      snowevap(j,lat) = snowevap(j,lat) 
     &      + snowev(j)*deltim*(hvap+hfus)/hvap
#else
        snowevap(j,lat) = snowevap(j,lat) + snowev(j)*deltim
#endif
#endif
#ifdef NFDR
#ifdef ICE
        es = fpvs(t2m(j,lat))
#else
        es = fpvs0(t2m(j,lat))
#endif
        qs = eps * es / (psexp(j) + epsm1*es)
        rh2m(j) = max (0.0e0, min (1.0e0, q2m(j,lat)/qs))*100.
        rhmax(j,lat) = max(rhmax(j,lat),rh2m(j))
        rhmin(j,lat) = min(rhmax(j,lat),rh2m(j))
        tmax24(j,lat) = max(tmax24(j,lat),t2m(j,lat))
        tmin24(j,lat) = min(tmin24(j,lat),t2m(j,lat))
#endif
        gflux(j,lat)=gflux(j,lat)+gflx(j)*dtf
        tmpmax(j,lat) = max(tmpmax(j,lat),t2m(j,lat))
        tmpmin(j,lat) = min(tmpmin(j,lat),t2m(j,lat))
        alhtfl(j,lat) = alhtfl(j,lat) + evap(j) * dtf
#ifdef HYDRO
        evcnp(j,lat) = evcnp(j,lat) + evcn(j) * dtf
#else
#ifdef NOAHYDRO
        if (slmsk(j,lat).eq.1) 
     &   evcnp(j,lat) = evcnp(j,lat) + evcn(j) * deltim
#endif
#endif
        ep(j,lat) = ep(j,lat) + ep1d(j) * dtf
#ifndef OSULSM1
        uustar(j,lat) = ustar(j)
        ffmm(j,lat) = fm(j)
        ffhh(j,lat) = fh(j)
#endif
      enddo
c
c     compute coefficient of evaporation in evapc
c
      do j=1,LONS2
        tstar(j) = tsea(j,lat)
        if (evapc(j).gt.1. e 0) evapc(j)=1.0 e 0
      enddo
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(tstar,LONS2,
     1            1,1,1,'tstar in gbphys')
      call maxmin(qss,LONS2,
     1            1,1,1,'qss in gbphys')
      call maxmin(evcn,LONS2,
     1            1,1,1,'evcn in gbphys')
#ifdef MP
      endif
#endif
#endif
c
c     over snow cover or ice or sea, coef of evap =1.0. e 0
c
      do 330 j=1,LONS2
      if ((sheleg(j,lat).gt.0. e 0) .or.
     1     (slmsk(j,lat).ne.1. e 0))
     2      evapc(j)=1. e 0
  330 continue
c
c  do vertical diffusion
c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            dtdt,psexp,-dtf,nmtvrdf,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            dtdt,-dtf,kdtvrdf,gda)
#endif
#ifdef VIC
      do i = 1, LONF2S
         if (tsea(i,lat).lt.273.15) then
           ls = (677.-0.07*(tsea(i,lat)-273.15))*4.1868*1000.0
           evp(i) = evap(i)/(rho(i)*ls)
         else
           le = (2.501-0.002361*(tsea(i,lat)-273.15))*1.e6
           evp(i) = evap(i)/(rho(i)*le)
         endif
         heat(i) = hflx(i)/rho(i)/cp_
      enddo
#endif
#ifdef CLDT93
      do k=1,levs_
      do j=1,LONS2
        hdf(j,k)=dtdt(j,k)
        qdf(j,k)=dqdt(j,k)
      end do 
      end do 
#endif
      call moninp(LONS2,LONF2S,
     1          levs_,nvdiff,dvdt,dudt,
     1         dtdt,dqdt,
     1         ugrs(1,1),vgrs(1,1),
     1         tgrs(1,1),qgrs(1,1),
     2         psexp,rb,cd,
     4         cdq,fm,fh,tstar,qss,evapc,wind,
c     5         si,del,sl,slk,rcl,deltim,lat,kdt,thour,
     5         sihyb,delhyb,slhyb,slkhyb,rcl,deltim,lat,kdt,thour,
     6         dusfc1,dvsfc1,dtsfc1,dqsfc1,
     7         hpbl(1,lat),gamt,gamq,slmsk(1,lat),
#ifdef ISOTOPE
     $         revap,dqsfc1i1,dqsfc1i2,qssi,tsea(1,lat),
     $         rrsvi1(1,1),rrsvi2(1,1),
#ifdef RSM      
     $         rlon(1,lat),rlat(1,lat),
#else
     $         xlon(1,lat),xlat(1,lat),
#endif
#endif
#ifdef CO2PROG
     $         co2fl(1,lat),dqsfc1i1,
#endif
#ifndef VIC
     $         kpbl)
#else
     $         kpbl,heat,evp)
#endif
#ifdef CLDT93
      do k=1,levs_
      do j=1,LONS2
        hdf(j,k)=dtdt(j,k)-hdf(j,k) 
        qdf(j,k)=dqdt(j,k)-qdf(j,k)
      end do 
      end do 
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,' monin'
      call maxmin(dudt(1,level),LONS2,
     1            1,1,numlevs,'dudt')
      call maxmin(dvdt(1,level),LONS2,
     1            1,1,numlevs,'dvdt')
      call maxmin(dtdt(1,level),LONS2,
     1            1,1,numlevs,'dtdt')
      call maxmin(dqdt(1,level),LONS2,
     1            1,1,numlevs,'dqdt')
#ifdef ISOTOPE
      call maxmin(dqdt(1,level+2*levs_),LONS2,
     1            1,1,numlevs,'dqdt2')
      call maxmin(dqdt(1,level+3*levs_),LONS2,
     1            1,1,numlevs,'dqdt3')
#endif
      call maxmin(dqsfc1(1),LONS2,
     1            1,1,1,'dqsfc1')
      call maxmin(dtsfc1(1),LONS2,
     1            1,1,1,'dtsfc1')
#ifdef MP
      endif
#endif
#endif
#ifdef DBG
#ifndef MP
      if(lat.eq.LATG2S)
     1 call rmaxmin(dqdt(1,1),LONS2,
     1              ktotal,1,
     1              ktotal,'dqdt after monin')
#endif
#endif
c
      do 340 j=1,LONS2
      dusfc(j,lat)=dusfc(j,lat)+dusfc1(j)*dtf
      dvsfc(j,lat)=dvsfc(j,lat)+dvsfc1(j)*dtf
      dtsfc(j,lat)=dtsfc(j,lat)+dtsfc1(j)*dtf
      dqsfc(j,lat)=dqsfc(j,lat)+dqsfc1(j)*dtf
#ifdef ISOTOPE      
      dqsfci1(j,lat)=dqsfci1(j,lat)+dqsfc1i1(j)*dtf
      dqsfci2(j,lat)=dqsfci2(j,lat)+dqsfc1i2(j)*dtf
#endif
#ifdef CO2PROG
      dqsfci1(j,lat)=dqsfci1(j,lat)+dqsfc1i1(j)*dtf/2500000 
#endif
340   continue
c
c  gravity wave drag
c
      do j = 1,LONS2
        var(j)=hprime(j,lat,1)
      enddo
#ifndef KIMGWD
      call gwdps(LONS2,LONF2S,
     1 levs_,dvdt,dudt,
     1 ugrs(1,1),vgrs(1,1),
     2 tgrs(1,1),qgrs(1,1),
     3  psexp,sihyb,delhyb,clhyb,slhyb,rcl,deltim,lat,kdt,
     $  hprime(1,lat,1),
     4  dusfcg,dvsfcg)
#else
      do j = 1, LONLS2
        oc(j) = hprime(j,lat,2)
      enddo
      do k = 1, 4
        do j = 1,LONLS2
          oa4(j,k) = hprime(j,lat,k+2)
          ol4(j,k) = hprime(j,lat,k+6)
        enddo
      enddo
      call gwd99(LONS2,LONF2S,
     1       levs_,dvdt,dudt,
     1       ugrs(1,1),vgrs(1,1),
     2       tgrs(1,1),qgrs(1,1),
     3       psexp,si,del,sl,slk,rcl,deltim,lat,kdt,
     4       var,oc,oa4,ol4,dusfcg,dvsfcg,kpbl,dxmeter)
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,' gwd'
      call maxmin(dudt(1,level),LONS2,1,1,numlevs,'dudt')
      call maxmin(dvdt(1,level),LONS2,1,1,numlevs,'dvdt')
#ifdef MP
      endif
#endif
#endif
c
      do 341 j=1,LONS2
      dugwd(j,lat)=dugwd(j,lat)+dusfcg(j)*dtf
      dvgwd(j,lat)=dvgwd(j,lat)+dvsfcg(j)*dtf
341   continue
c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            dudt,psexp, dtf,nmuvrdf,lat)
      call znlacm(LONS2,LONF2S,
     1            dvdt,psexp, dtf,nmvvrdf,lat)
      call znlacm(LONS2,LONF2S,
     1            dtdt,psexp, dtf,nmtvrdf,lat)
      call znlacm(LONS2,LONF2S,
     1            dqdt,psexp, dtf,nmqvrdf,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            dudt, dtf,kduvrdf,gda)
      call accdia(LONS2,LONF2S,
     1            dvdt, dtf,kdvvrdf,gda)
      call accdia(LONS2,LONF2S,
     1            dtdt, dtf,kdtvrdf,gda)
      call accdia(LONS2,LONF2S,
     1            dqdt, dtf,kdqvrdf,gda)
#endif
#ifndef SCM
      if(inistp.le.1) then
c
c   get vertical motion (cb/sec) in vvel
c
#ifdef RSMVECTORIZE
#define OMEGAS omegasv
#else
#define OMEGAS omegas
#endif
      call OMEGAS(LONS2,LONF2S,
     1            levs_,
     1            pphigr(1),plamgr(1),w2,
     1            ugrs(1,1),vgrs(1,1),
#ifndef RSM
     2            xgrs(1,1),delhyb,rbs2(lat),vvel,
#else
     2            xgrs(1,1),delhyb,rbs2,vvel,
#endif
     3            psexp,slhyb)
      endif
#else
      print *, ' OMEGA : vvel -> xgrs in SCM '
#endif
c
c   w2 is used to store wind speed
c
      do k=1,levs_
        do j=1,LONS2
          w2(j,k)=rcs*sqrt(ugrs(j,k)**2+vgrs(j,k)**2)
        enddo
      enddo
c
#ifdef DBG
#ifndef MP
      if(lat.eq.LATG2S)
     1 call rmaxmin(qgrs(1,1),LONS2,
     1              ktotal,1,ktotal,'gq0 before monin')
#endif
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
#ifdef ISOTOPE
      call maxmin(qgrs(1,level),LONS2,
     1            1,1,numlevs,'gq0 before monin')
      call maxmin(qgrs(1,level+2*levs_),LONS2,
     1            1,1,numlevs,'gq2 before monin')
      call maxmin(qgrs(1,level+3*levs_),LONS2,
     1            1,1,numlevs,'gq3 before monin')
#endif
#ifdef MP
      endif
#endif
#endif
      do k=1,levs_
        do j=1,LONS2
          gt0(j,k)=tgrs(j,k)+dtdt(j,k)*dt2
          gu0(j,k)=ugrs(j,k)+dudt(j,k)*dt2
          gv0(j,k)=vgrs(j,k)+dvdt(j,k)*dt2
        enddo
      enddo
      do k=1,levh_
        do j=1,LONS2
          gq0(j,k)=qgrs(j,k)+dqdt(j,k)*dt2
        enddo
      enddo
#ifdef DBG
#ifndef MP
       if(lat.eq.LATG2S)
     1 call rmaxmin(gq0(1,1),LONS2,
     1              ktotal,1,ktotal,'gq0 after monin')
#endif
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
#ifdef ISOTOPE
      call maxmin(gq0(1,level),LONS2,
     1            1,1,numlevs,'gq0 after monin')
      call maxmin(gq0(1,level+2*levs_),LONS2,
     1            1,1,numlevs,'gq2 after monin')
      call maxmin(gq0(1,level+3*levs_),LONS2,
     1            1,1,numlevs,'gq3 after monin')
#endif
#ifdef MP
      endif
#endif
#endif
c
#ifndef O3C
      if(ntrace .ge. 2) then
        call o3phys(LONS2,LONF2S,
     1              levs_,deltim,pgr,gq0(1,ktrace),
#ifndef RSM
     1              gq0(1,ktrace),xlat(1,lat),
#else
     1              gq0(1,ktrace),rlat(1,lat),
#endif
     2              slhyb)
      endif
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(gt0(1,level),LONS2,
     1            1,1,numlevs,'gt0')
      call maxmin(gq0(1,level),LONS2,
     1            1,1,numlevs,'gq0')
      call maxmin(gu0(1,level),LONS2,
     1            1,1,numlevs,'gu0')
      call maxmin(gv0(1,level),LONS2,
     1            1,1,numlevs,'gv0')
#ifdef ISOTOPE
      call maxmin(gq0(1,level+2*levs_),LONS2,
     1            1,1,numlevs,'gq2')
      call maxmin(gq0(1,level+3*levs_),LONS2,
     1            1,1,numlevs,'gq3')
#endif
#ifdef MP
      endif
#endif
#endif
c
c  call to gwater
c
      if(inistp.le.1) then
c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,gt0(1,1),
     1            psexp,-frain,nmtconv,lat)
      call znlacm(LONS2,LONF2S,gq0(1,1),
     1            psexp,-frain,nmqconv,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1),-frain,kdtconv,gda)
      call accdia(LONS2,LONF2S,
     1            gq0(1,1),-frain,kdqconv,gda)
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
c
c store q & ql integrated vertically and rain
c
      qck1(:,:)=0.0  ! q, qc, qw, 
      rck1(:)=0.0    ! rain
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck1(i,1)=qck1(i,1) + gq0(i,k)*ds   ! q
        end do
      end do
      if (ncloud.gt.0) then
      do n=1,ncloud
      do k=1,levs_
        kk=kcloud+(n-1)*levs_+k-1
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck1(i,n+1)=qck1(i,n+1) + gq0(i,kk)*ds   ! qc,qw,
        end do
      end do
      end do
      endif
      do n=1,ncloud+1
      do i=1,LONS2
        qck1(i,n)=qck1(i,n) * psexp(i) * ucnv1   ! [kg/kg] -> [kg/m2]=[mm]
      end do
      end do

#ifdef MP
      endif
#endif
#endif

#ifdef SAS
#ifdef SCM
      call sascnv(LONS2,LONF2S,levs_,jcapset,deltim,del,sl,slk,psexp,
     2            gq0(1,1),gt0(1,1),cld1d,
     3            rain1,kbot,ktop,kuo,w2,lat,slmsk(1,lat),xgrs,
#else
      call sascnv(LONS2,LONF2S,
#ifdef RSM
     1            levs_,jcapr,deltim,delhyb,slhyb,slkhyb,psexp,
#else
     1            levs_,jcap_,deltim,delhyb,slhyb,slkhyb,psexp,
#endif
     2            gq0(1,1),gt0(1,1),cld1d,
     3            rain1,kbot,ktop,kuo,w2,lat,slmsk(1,lat),vvel,
#endif
#ifdef UV_DIFF
     &            gu0(1,1),gv0(1,1),rcs,
#endif
#ifdef UV_TRAN
     &            gu0(1,1),gv0(1,1),rcs,
#endif
#ifdef MUL_CLDTOP
     &            xkt2,NCLDTOP,
#endif
#ifdef CLM_CWF
     &            cgs2,stepone,lwconv,
#endif
#ifdef SAS_DIAG
     +            dcu,dcv,dct,dcq,dch,fcu,fcd,
     +            deltb,delqb,delhb,cbmf,
#endif
     4            hpbl(1,lat),gamt,gamq,qcicnv,qrscnv,dxmeter)
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'sascnv'
      call maxmin(gt0(1,level),LONS2,
     1            1,1,numlevs,'gt0')
      call maxmin(gq0(1,level),LONS2,
     1            1,1,numlevs,'gq0')
      call maxmin(gu0(1,level),LONS2,
     1            1,1,numlevs,'gu0')
      call maxmin(gv0(1,level),LONS2,
     1            1,1,numlevs,'gv0')
#ifdef MP
      endif
#endif
#endif
#endif
#ifdef RAS
      dtras = deltim / frain
      il1r = 1 + 8 * LONF2S
      il2r = il1r + 10 * LONF2S * lm
      il3r = il2r +  LONF2S * (5*lm+16)
      clw(1:LONS2,1:levs_)=0.0
      cmf(1:LONS2,1:levs_)=0.0
      detr(1:LONS2,1:levs_)=0.0
#ifdef ISOTOPE
      do j=1,LONS2
        eqfo(j)=meqfo(j,lat)
        eqfd(j)=meqfd(j,lat)
      enddo
#endif
ckei ras is moved into longitudinal loop.
c      call rascnv(LONS2, LONF2S,
c     1            levs_, lm, nstrp, dtras, si, sl
c     *,           krmin, krmax, ncrnd, afac, rannum, ufac
c     *,           rgas, cp, grav, alhl
c     *,           gt0(1,1), gq0(1,1)
c     *,           gu0(1,1), gv0(1,1), psexp
c     *,           rain1, kbot, ktop, kuo, lat, cd
c     *,           sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns
c     *,           clw ,cmf,detr)
       do j=1,LONS2
c
         do k=1,levs_
           si0(k)=sihyb(j,k)
           sl0(k)=slhyb(j,k)
           del0(k)=delhyb(j,k)
         enddo
         si0(levp1_)=sihyb(j,levp1_)
         call setrasv1_1d(levs_,si0,sl0,del0,cp_,rd_
     *,             sig,sgb,prh,prj,hpk,fpk,ods,prns
     *,             lm,nstrp,rasal,krmin,krmax,afac,ufac)
c
         do k=1,levs_
           gt00(k)=gt0(j,k)
           gq00(k)=gq0(j,k)
           gu00(k)=gu0(j,k)
           gv00(k)=gv0(j,k)
           clw0(k)=0.
           cmf0(k)=0.
           detr0(k)=0.
         enddo
#ifdef ISOTOPE         
         do k=levs_+1,levh_
           gq00(k)=gq0(j,k)
         enddo
#endif
c
c#ifdef ISOTOPE         
c         do k=1,2
c            tmpr(k)=rain1i(j,k)
c            tmpc(k)=cond1i(j,k)
c         enddo
c#endif
         call rascnp_1d(1,1,
     1            levs_, lm, nstrp, dtras, si0, sl0
     *,           krmin, krmax, ncrnd, afac, rannum, ufac
     *,           rgas, cp, grav, alhl
     *,           gt00, gq00, gu00, gv00, psexp(j)
     *,           rain1(j), kbot(j), ktop(j), kuo(j), lat, cd(j)
     *,           sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns
#ifdef ISOTOPE
     *,           ntrac_,tmpr(1),tmpc(1),cond1(j)
     *,           eqfo(j),eqfd(j)
#endif
     *,           clw0 ,cmf0,detr0)
#ifdef ISOTOPE         
         do k=1,2
            rain1i(j,k)=tmpr(k)
            cond1i(j,k)=tmpc(k)
         enddo
#endif
         do k=1,levs_
           gt0(j,k)=gt00(k)
           gq0(j,k)=gq00(k)
           gu0(j,k)=gu00(k)
           gv0(j,k)=gv00(k)
           clw(j,k)=clw0(k)
           cmf(j,k)=cmf0(k)
           detr(j,k)=detr0(k)
         enddo
#ifdef ISOTOPE
         do k=levs_+1,levh_
           gq0(j,k)=gq00(k)
         enddo
#endif
       enddo
c
       do j=1,LONS2
         cld1d(j)=0.
       enddo
#endif
#ifdef RASV2
       dtras = deltim / frain
       do j=1,LONS2
         do k=1,levs_
           prsi(k)=sihyb(j,k)*psexp(j)
           prsl(k)=slhyb(j,k)*psexp(j)
         enddo
         prsi(levs_+1)=sihyb(j,levs_+1)*psexp(j)
c
         call setrasv2_hyb(levs_,prsi,prsl,cp_,rd_
     *,             sig,sgb,prj
     *,             krmin,krmax,kfmax
     *,             mct,kctop)
c
         do k=1,levs_
           gt00(k)=gt0(j,k)
           gq00(k)=gq0(j,k)
           gu00(k)=gu0(j,k)
           gv00(k)=gv0(j,k)
         enddo
         call rascnv2(levs_, dtras
     *,               krmin, krmax, kfmax, ncrnd, rannum, pdd
     *,               mct, kctop
     *,               rgas, cp, grav, alhl, hfus
     *,               gt00, gq00, gu00, gv00, psexp(j), kpbl(j)
     *,               rain1(j), kbot(j), ktop(j), kuo(j), lat, cd(j)
     *,               sig, prj, sgb, rasal, clw0, clt, dsfc(j))
         do k=1,levs_
           gt0(j,k)=gt00(k)
           gq0(j,k)=gq00(k)
           gu0(j,k)=gu00(k)
           gv0(j,k)=gv00(k)
           clw(j,k)=clw0(levs_-k+1)
         enddo
       enddo
       do j=1,LONS2
         cld1d(j)=0.
       enddo
#endif
#ifdef CCMCNV
      call ccmcnv_driver(LONS2,LONF2S,
     1            levs_,
#ifdef SCM
     1            jcapset,deltim,si,del,sl,slk,psexp,kdt,
     2            gq0(1,1),gt0(1,1),cld1d,
     3            rain1,kbot,ktop,kuo,w2,lat,slmsk(1,lat),xgrs,
#else
     1            jcap_,deltim,sihyb,delhyb,slhyb,slkhyb,psexp,kdt,
     2            gq0(1,1),gt0(1,1),cld1d,
     3            rain1,kbot,ktop,kuo,w2,lat,slmsk(1,lat),vvel,
#endif
     4            hpbl(1,lat),gamt,gamq,qcicnv,qrscnv,dxmeter)
#endif

#ifdef ZM2
      if (shour.eq.0.0.or.idpinizm(lat).eq.0) then
        do k = 1, levs_
          do i = 1, LONS2
            tpzm(i,k,lat) = gt0(i,k)
            qpzm(i,k,lat) = max(gq0(i,k),epsq)
            tp1zm(i,k,lat) = gt0(i,k)
            qp1zm(i,k,lat) = max(gq0(i,k),epsq)
          enddo
        enddo
        idpinizm(lat)=1
      endif

      call zm_conv_driver(LONS2,LONF2S,
     1            levs_,
#ifdef SCM
     1            jcapset,
#else
     1            jcap_,
#endif
     1            deltim,sihyb,psexp,kdt,
     2            gq0(1,1),gt0(1,1),
     3            rain1,kbot,ktop,kuo,lat,
     4            hpbl(1,lat),gamt,gamq,qcicnv,qrscnv,dxmeter,
     5            tpzm(1,1,lat),qpzm(1,1,lat),
     6            tp1zm(1,1,lat),qp1zm(1,1,lat),clw,detr)
       cld1d(1:LONS2)=0.
#endif

#ifdef KF2
      call kf2cnv(LONS2,LONF2S,levs_,gu0(1,1),gv0(1,1),
     1           gt0(1,1),gq0(1,1),sihyb,slhyb,delhyb,psexp,
c###
     2           ncloud,clw,cli,
c###
#ifdef SCM
     2           xgrs,deltim,dxmeter,rain1,kbot,ktop,kuo,lat)
#else
     2           vvel,deltim,dxmeter,rain1,kbot,ktop,kuo,lat)
#endif
      do j=1,LONS2
        cld1d(j)=0.
      enddo
#endif
c
      do  j=1,LONS2
         cldwrk(j,lat) = cldwrk(j,lat) + cld1d(j) * dtf
      enddo

c
c!!! treatment of detrained cloud water/ice
#ifdef CLD1
      ndc=0
#endif
#ifdef CLD2
      ndc=1
#endif
#ifdef CLDT93
      ndc=0
#endif
#ifdef CLDWSMMPS
      ndc=1
#endif
#ifdef CLD3
      ndc=1
#endif
#ifdef CLD5
      ndc=2
#endif
#ifdef CLD6
      ndc=2
#endif
      if (ndc.eq.1) then  ! add cloud water/ice
#ifdef RAS
        do k = 1, levs_
          do i = 1, LONS2
            gq0(i,kcloud+k-1)=gq0(i,kcloud+k-1)+clw(i,k)
          end do     
        end do     
#endif
#ifdef KF2
        do k = 1, levs_
          do i = 1, LONS2
            gq0(i,kcloud+k-1)=gq0(i,kcloud+k-1)+clw(i,k)
          end do     
        end do     
#endif
#ifdef ZM2
        do k = 1, levs_
          do i = 1, LONS2
            gq0(i,kcloud+k-1)=gq0(i,kcloud+k-1)+clw(i,k)
          end do
        end do
#endif
      else if (ndc.eq.2) then ! add cloud water,ice
#ifdef RAS
c     !!! divide clw into water and ice like cldprpi.F
        do k = 1, levs_
          do i = 1, LONS2
            dtemp = 263.16 - gt0(i,k)
            fice = amin1(1.0, amax1(0.0e0, dtemp*0.05e0 ))
            clw1 = clw(i,k) * fice  ! ice
            clw2 = clw(i,k) - clw1  ! water
            gq0(i,kcloud+k-1)      = gq0(i,kcloud+k-1)+clw2
            gq0(i,kcloud+levs_+k-1)= gq0(i,kcloud+levs_+k-1)+clw1
          end do     
        end do     
#endif
#ifdef ZM2
c     !!! divide clw into water and ice like cldprpi.F
        do k = 1, levs_
          do i = 1, LONS2
            dtemp = 263.16 - gt0(i,k)
            fice = amin1(1.0, amax1(0.0e0, dtemp*0.05e0 ))
            clw1 = clw(i,k) * fice  ! ice
            clw2 = clw(i,k) - clw1  ! water
            gq0(i,kcloud+k-1)      = gq0(i,kcloud+k-1)+clw2
            gq0(i,kcloud+levs_+k-1)= gq0(i,kcloud+levs_+k-1)+clw1
          end do     
        end do     
#endif
#ifdef KF2
        do k = 1, levs_
          do i = 1, LONS2
            gq0(i,kcloud+k-1)=gq0(i,kcloud+k-1)+clw(i,k)              ! water
            gq0(i,kcloud+levs_+k-1)= gq0(i,kcloud+levs_+k-1)+cli(i,k) ! ice
          end do     
        end do     
#endif
      endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
c
c store q & ql integrated vertically and rain
c
      qck2(:,:)=0.0  ! q, qc, qw, 
      rck2(:)=0.0    ! rain
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck2(i,1)=qck2(i,1) + gq0(i,k)*ds   ! q
        end do
      end do
      if (ncloud.gt.0) then
      do n=1,ncloud
      do k=1,levs_
        kk=kcloud+(n-1)*levs_+k-1
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck2(i,n+1)=qck2(i,n+1) + gq0(i,kk)*ds   ! qc,qw,
        end do
      end do
      end do
#ifdef CLDT93
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck2(i,2)=qck2(i,2) + clw(i,k)*ds   ! add detrained cloud water
        end do     
      end do     
#endif
      endif
      do n=1,ncloud+1
      do i=1,LONS2
        qck2(i,n)=qck2(i,n) * psexp(i) * ucnv1   !([kg/kg])-> [kg/m2]=[mm]
      end do
      end do
      do i=1,LONS2
        rck2(i)=rain1(i)               * ucnv2   ! [m]     -> [kg/m2]=[mm]
      end do
c
c compare
c
      do i=1,LONS2
        qsum1=0.0
        qsum2=0.0
        do n=1,ncloud+1
          qsum1=qsum1+qck1(i,n)
          qsum2=qsum2+qck2(i,n)
        end do
        qsum1=qsum1+rck1(i)
        qsum2=qsum2+rck2(i)
        if (qsum1.ne.qsum2) then
          qsum3=(qsum2-qsum1)/qsum1*100.0  ! %
          if (abs(qsum3).gt.1.0) then
            print '(a7,x,i3,x,f5.2,2(x,f10.6))'
     *            ,'qck-1:',i,qsum3,qsum1,qsum2
            do n=1,ncloud+1
              print '(a7,x,i2,3(x,f10.6))'
     *            ,'qck2-1',n,qck2(i,n)-qck1(i,n),qck1(i,n),qck2(i,n)
            end do
            print '(a10,x,f10.6)'
     *            ,'qck2-1(r)',rck2(i)
          endif
        endif
      end do
#ifdef MP
      endif
#endif
#endif

c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            gt0(1,1),psexp, frain,nmtconv,lat)
      call znlacm(LONS2,LONF2S,
     1            gq0(1,1),psexp, frain,nmqconv,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1), frain,kdtconv,gda)
      call accdia(LONS2,LONF2S,
     1            gq0(1,1), frain,kdqconv,gda)
#endif
c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            gt0(1,1),psexp,-frain,nmtshal,lat)
      call znlacm(LONS2,LONF2S,
     1            gq0(1,1),psexp,-frain,nmqshal,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1),-frain,kdtshal,gda)
      call accdia(LONS2,LONF2S,
     1            gq0(1,1),-frain,kdqshal,gda)
#endif
#ifdef CCMSCV
      call ccmscv_driver(LONS2,LONF2S,
     1            levs_,
#ifdef SCM
     1            jcapset,deltim,sihyb,delhyb,slhyb,slkhyb,psexp,kdt,
     2            gq0(1,1),gt0(1,1),
     3            rain2,kbot,ktop,kuo,w2,lat,slmsk(1,lat),xgrs,
#else
     1            jcap_,deltim,sihyb,delhyb,slhyb,slkhyb,psexp,kdt,
     2            gq0(1,1),gt0(1,1),
     3            rain2,kbot,ktop,kuo,w2,lat,slmsk(1,lat),vvel,
#endif
     4            hpbl(1,lat),gamt,gamq,qcicnv,qrscnv,dxmeter)
      do j=1,LONS2
        rainc(j)=frain*(rain1(j)+rain2(j))
        bengsh(j,lat)=bengsh(j,lat)+rainc(j)
      enddo
#else
      call shalcv(LONS2,LONF2S,
     1            levs_,levshc,nvdiff,deltim,
     $            kuo,psexp,gq0(1,1),gt0(1,1),
     1            delhyb,slhyb,slkhyb,sihyb)
      do j=1,LONS2
        rainc(j)=frain*rain1(j)
        bengsh(j,lat)=bengsh(j,lat)+rainc(j)
#ifdef ISOTOPE
        condc(j)=frain*cond1(j)
        bengshc(j,lat)=bengshc(j,lat)+condc(j)
        meqfo(j,lat)=min(meqfo(j,lat),eqfo(j))
        meqfd(j,lat)=min(meqfd(j,lat),eqfd(j))
#endif
      enddo
#ifdef ISOTOPE
      do j=1,LONS2
        do k=1,2
          rainci(j,k)=frain*rain1i(j,k)
          condci(j,k)=frain*cond1i(j,k)
          if (k.eq.1) then
            bengshi1(j,lat)=bengshi1(j,lat)+rainci(j,k) !! [m]
            bengshci1(j,lat)=bengshci1(j,lat)+condci(j,k) !! [m]
          elseif (k.eq.2) then
            bengshi2(j,lat)=bengshi2(j,lat)+rainci(j,k)
            bengshci2(j,lat)=bengshci2(j,lat)+condci(j,k)
          endif
        enddo
      enddo
#endif
#endif
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'shalcv'
      call maxmin(gt0(1,level),LONS2,
     1            1,1,numlevs,'gt0')
      call maxmin(gq0(1,level),LONS2,
     1            1,1,numlevs,'gq0')
      if(ntrace.ge.2) then
      call maxmin(gq0(1,ktrace),LONS2,
     1            1,1,levs_,'trace')
      endif
      if(ncloud.ge.1) then
      call maxmin(gq0(1,kcloud),LONS2,
     1            1,1,levs_,'cloud')
      endif
#ifdef MP
      endif
#endif
#endif
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            gt0(1,1),psexp, frain,nmtshal,lat)
      call znlacm(LONS2,LONF2S,
     1            gq0(1,1),psexp, frain,nmqshal,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1), frain,kdtshal,gda)
      call accdia(LONS2,LONF2S,
     1            gq0(1,1), frain,kdqshal,gda)
#endif
c
#ifndef CLDSLINGO
#define CNVCLD cnvcld
#else
#define CNVCLD cnvcldn
#endif
      call CNVCLD(LONS2,LONF2S,
     1            clstp,rainc,kbot,ktop,
#ifndef RSM
     1            cv(1,lat),cvb(1,lat),cvt(1,lat),count)
#else
     1            cv(1,lat),cvb(1,lat),cvt(1,lat),1.)
#endif
c
#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            gt0(1,1),psexp,-frain,nmtlarg,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1),-frain,kdtlarg,gda)
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
c
c store q & ql integrated vertically and rain
c
      qck1(:,:)=0.0  ! q, qc, qw, 
      rck1(:)=0.0    ! rain
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck1(i,1)=qck1(i,1) + gq0(i,k)*ds   ! q
        end do
      end do
      if (ncloud.gt.0) then
      do n=1,ncloud
      do k=1,levs_
        kk=kcloud+(n-1)*levs_+k-1
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck1(i,n+1)=qck1(i,n+1) + gq0(i,kk)*ds   ! qc,qw,
        end do
      end do
      end do
#ifdef CLDT93
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck1(i,2)=qck1(i,2) + clw(i,k)*ds   ! add detrained cloud water
        end do     
      end do     
#endif
      endif
      do n=1,ncloud+1
      do i=1,LONS2
        qck1(i,n)=qck1(i,n) * psexp(i) * ucnv1   ! [kg/kg] -> [kg/m2]=[mm]
      end do
      end do

#ifdef MP
      endif
#endif
#endif

#ifdef CLD1
      call lrgscl(LONS2,LONF2S,
     1            levs_,deltim,psexp,
     1            gt0(1,1),gq0(1,1),slhyb,delhyb,slkhyb,
#ifdef CONDENS
     +             qcilrg,qrslrg,
#endif
#ifdef ISOTOPE
     &            ntrac_,rain1i,cond1i,cond1,
#endif
#ifndef RASV2
     &              rain1,lat)
#else
     &              rain1,lat,clw)
#endif
#endif

#ifdef CLD2
c
c  check for first time step of the segment run
c
cc      if (shour.eq.0.) then
      if (shour.eq.0.0.or.idpini(lat).eq.0) then
        do k = 1, levs_
          do i = 1, LONS2
            tp(i,k,lat) = gt0(i,k)
            qp(i,k,lat) = max(gq0(i,k),epsq)
            tp1(i,k,lat) = gt0(i,k)
            qp1(i,k,lat) = max(gq0(i,k),epsq)
          enddo
        enddo
        do i = 1, LONS2
          psp(i,lat)  = psexp(i)
          psp1(i,lat) = psexp(i)
        enddo
        idpini(lat)=1
      endif
      call gscond(LONS2,LONF2S,
     1            levs_,deltim,delhyb,slhyb,psexp,
     &            gq0(1,1),gq0(1,kcloud),gt0(1,1),lat,
     &            tp(1,1,lat),qp(1,1,lat),psp(1,lat),
     &            tp1(1,1,lat),qp1(1,1,lat),psp1(1,lat),rcs)
      call precpd(LONS2,LONF2S,
     &            levs_,deltim,delhyb,slhyb,psexp,
     &            gq0(1,1),gq0(1,kcloud),gt0(1,1),rain1,lat,rcs)
#endif

#ifdef CLDWSMMPS
      call cloud2(LONS2,LONF2S,
     1            levs_,deltim,psexp,
     1            gt0(1,1),gq0(1,1),gq0(1,kcloud),ncloud,
     $            slhyb,delhyb,slkhyb,
#ifdef SCM
     2            rain1,lat,xgrs,kdt,fhour)
#else
     2            rain1,lat,vvel,kdt,fhour)
#endif
#endif

#ifdef CLD3
      call cloud3(LONS2,LONF2S,
     1            levs_,deltim,psexp,
     1            gt0(1,1),gq0(1,1),gq0(1,kcloud),ncloud,
     $            slhyb,delhyb,slkhyb,
#ifdef SCM
     2            rain1,lat,xgrs,kdt,fhour)
#else
     2            rain1,lat,vvel,kdt,fhour)
#endif
#endif

#ifdef CONDENS
      do k = 1, levs_
      do i = 1, LONF2S
        qcilrg(i,k) = gq0(i,kcloud+k-1)
        qrslrg(i,k) = gq0(i,kcloud+levs_+k-1)
      enddo
      enddo
#endif

#ifdef CLD5
      call cloud5(LONS2,LONF2S,
     1            levs_,deltim,psexp,
     1            gt0(1,1),gq0(1,1),gq0(1,kcloud),ncloud,
     $            slhyb,delhyb,slkhyb,
#ifdef SCM
     2            rain1,lat,xgrs,kdt,fhour)
#else
     2            rain1,lat,vvel,kdt,fhour)
#endif
#endif

#ifdef CLD6
      call cloud6(LONS2,LONF2S,
     1            levs_,deltim,psexp,
     1            gt0(1,1),gq0(1,1),gq0(1,kcloud),ncloud,
     $            slhyb,delhyb,slkhyb,
#ifdef SCM
     2            rain1,lat,xgrs,kdt,fhour)
#else
     2            rain1,lat,vvel,kdt,fhour)
#endif
#endif

#ifdef CLDT93
          kcamn = 2*levs_+1  ! suppose cloud amount is 3rd in tracer
          call t93clds_driver(LONS2,LONF2S,
     1                  levs_,deltim,sihyb,slhyb,delhyb,psexp,
     2                  gq0(1,1),gq0(1,kcloud),gt0(1,1),gq0(1,kcamn),
     3                  hdf,qdf,dqsfc1,dtsfc1,hpbl(1,lat),slmsk(1,lat),
#ifdef SCM
     4                  xgrs,clw,hsw,hlw(1,1,lat),
#else
     4                  vvel,clw,hsw,hlw(1,1,lat),
#endif
     5                  cmf,detr,solin(1,1,lat),
     &                  rain1,lat)
#endif

#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
c
c store q & ql integrated vertically and rain
c
      qck2(:,:)=0.0  ! q, qc, qw, 
      rck2(:)=0.0    ! rain
      do k=1,levs_
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck2(i,1)=qck2(i,1) + gq0(i,k)*ds   ! q
        end do
      end do
      if (ncloud.gt.0) then
      do n=1,ncloud
      do k=1,levs_
        kk=kcloud+(n-1)*levs_+k-1
        do i=1,LONS2
          ds=sihyb(i,k)-sihyb(i,k+1)
          qck2(i,n+1)=qck2(i,n+1) + gq0(i,kk)*ds   ! qc,qw,
        end do
      end do
      end do
      endif
      do n=1,ncloud+1
      do i=1,LONS2
        qck2(i,n)=qck2(i,n) * psexp(i) * ucnv1   !([kg/kg])-> [kg/m2]=[mm]
      end do
      end do
      do i=1,LONS2
        rck2(i)=rain1(i)               * ucnv2   ! [m]     -> [kg/m2]=[mm]
      end do
c
c compare
c
      do i=1,LONS2
        qsum1=0.0
        qsum2=0.0
        do n=1,ncloud+1
          qsum1=qsum1+qck1(i,n)
          qsum2=qsum2+qck2(i,n)
        end do
        qsum1=qsum1+rck1(i)
        qsum2=qsum2+rck2(i)
        if (qsum1.ne.qsum2) then
          qsum3=(qsum2-qsum1)/qsum1*100.0  ! %
          if (abs(qsum3).gt.1.0) then
            print '(a7,x,i3,x,f5.2,2(x,f10.6))'
     *            ,'qck-2:',i,qsum3,qsum1,qsum2
            do n=1,ncloud+1
              print '(a7,x,i2,3(x,f10.6))'
     *            ,'qck2-2',n,qck2(i,n)-qck1(i,n),qck1(i,n),qck2(i,n)
            end do
            print '(a10,x,f10.6)'
     *            ,'qck2-2(r)',rck2(i)
          endif
        endif
      end do
#ifdef MP
      endif
#endif
#endif

#ifdef DGZ
      call znlacm(LONS2,LONF2S,
     1            gt0(1,1),psexp, frain,nmtlarg,lat)
#endif
#ifdef DG3
      call accdia(LONS2,LONF2S,
     1            gt0(1,1), frain,kdtlarg,gda)
#endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'lrgscl'
      call maxmin(gt0(1,level),LONS2,
     1            1,1,numlevs,'gt0')
      call maxmin(gq0(1,level),LONS2,
     1            1,1,numlevs,'gq0')
      call maxmin(gu0(1,level),LONS2,
     1            1,1,numlevs,'gu0')
      call maxmin(gv0(1,level),LONS2,
     1            1,1,numlevs,'gv0')
#ifdef MP
      endif
#endif
#ifndef MP
       if(lat.eq.latg2_)
     1 call rmaxmin(gq0(1,1),LONS2,
     1              ktotal,1,ktotal,'gq0 after cloud')
#endif
#endif
c
      do j=1,LONS2
        rainl(j)=frain*rain1(j)
        rain(j)=rainc(j)+rainl(j)
        geshem(j,lat)=geshem(j,lat)+rain(j)        
#ifdef ISOTOPE
        condl(j)=frain*cond1(j)
        cond(j)=condc(j)+condl(j)
        geshemc(j,lat)=geshemc(j,lat)+cond(j)
#endif
#ifdef NOAHYDRO
        geshem2(j,lat)=geshem2(j,lat)+prcp(j,lat)
#endif
#ifdef NFDR
c
c  precipitation duration and 24 hour precipitation
c
        if(rain(j).gt.0.) then
          prcpdur(j,lat)=prcpdur(j,lat)+dtf
        endif
        ppt24(j,lat)=ppt24(j,lat)+rain(j)
#endif
      enddo
c
#ifdef ISOTOPE
      do j=1,LONS2
        do k=1,2
          rainli(j,k)=frain*rain1i(j,k)
          raini(j,k)=rainci(j,k)+rainli(j,k)
          if (k.eq.1) then
            geshemi1(j,lat)=geshemi1(j,lat)+raini(j,k) !! [m]
          elseif (k.eq.2) then
            geshemi2(j,lat)=geshemi2(j,lat)+raini(j,k)
          endif
          condli(j,k)=frain*cond1i(j,k)
          condi(j,k)=condci(j,k)+condli(j,k)
          if (k.eq.1) then
            geshemci1(j,lat)=geshemci1(j,lat)+condi(j,k) !! [m]
          elseif (k.eq.2) then
            geshemci2(j,lat)=geshemci2(j,lat)+condi(j,k)
          endif
        enddo
      enddo
      do k=1,levs_
        do j=1,LONS2
c          if ((mype.eq.23).and.(j.eq.60)) then
c            print'(a6,3i4,3e15.7)','gbCHK3',j,k,lat,
c     $            gq0(j,k),gq0(j,k+3*levs_),
c     $            gq0(j,k+3*levs_)/gq0(j,k)
c          endif
c          if ((gq0(j,k).le.1.D-5).and.(gq0(j,k).gt.0.D0).and.
c     $           ((gq0(j,k+2*levs_)/gq0(j,k).gt.1.05D0).or.
c     $            (gq0(j,k+3*levs_)/gq0(j,k).gt.1.40D0))) then
c            print'(a6,4i4,3e15.7)','gbCHK2',mype,j,k,lat,
c     $            gq0(j,k),
c     $            gq0(j,k+2*levs_)/gq0(j,k),
c     $            gq0(j,k+2*levs_)/gq0(j,k)
c          endif
        enddo
      enddo
#endif
c
c  estimate t850 for rain-snow decision
c
      do j = 1, LONS2
        pk(j) = psexp(j) * sl(1)
        t850(j) = gt0(j,1)
      enddo
      do k = 1, levs_ - 1
        do j = 1, LONS2
          pkp(j) = psexp(j) * sl(k+1)
          if(pk(j) .gt. p850 .and. pkp(j) .le. p850) then
            t850(j) = gt0(j,k) - (pk(j) - p850) /
     &                (pk(j) - pkp(j)) * (gt0(j,k)
     &                - gt0(j,k+1))
          endif
          pk(j) = pkp(j)
        enddo
      enddo
c
c  factor=weighted mean tep.
c
      do j=1,LONS2
#ifdef MRGLSM
c
c  save rain and (t850 vs 273) to prcp and srflag; do not update sheleg
c
       prcp(j,lat) = rain(j)
       srflag(j,lat) = 0.
       if(t850(j).le.critsnow) srflag(j,lat)=1.
#else
c  update sheleg
        if(t850(j).le.critsnow) then
          snowfl(j)=rain(j)
          if(slmsk(j,lat).ne.0.) then
            sheleg(j,lat)=sheleg(j,lat)+1.e3*rain(j)
            snowfall(j,lat)=snowfall(j,lat)+rain(j)
          endif
          rain(j)=0.
        else
          snowfl(j)=0.
        endif
#endif
ckei
#ifdef RSM        
        latj = (lat-1)*LONF2S
        if (rltb(j+latj).ne.1) then
          bengsh(j,lat)=0.
          geshem(j,lat)=0.
#ifdef ISOTOPE
          bengshi1(j,lat)=0.
          geshemi1(j,lat)=0.
          bengshi2(j,lat)=0.
          geshemi2(j,lat)=0.
#endif
#ifdef NOAHYDRO
          geshem2(j,lat)=0.
#endif
          snowfall(j,lat)=0.
          prcp(j,lat)=0.
          srflag(j,lat)=0.
c          if (mype.eq.master) then
c            print*,'location',lat,j,rlon(j,lat),rlat(j,lat)
c          endif
        endif
#endif
      enddo

#ifndef MRGLSM
c
c  update soil moisture and canopy water after precipitation has been
c  computed
c
      call progt2(LONS2,LONF2S,
     1  lsoil_,
     &  rhscnpy,rhsmc,ai,bi,cci,smsoil,
     &  slmsk(1,lat),canopy(1,lat),rain,runof,snowmt,
#ifdef HYDRO
     &  zsoil,soiltyp,sigmaf,deltim,lat,hydrow)
#else
     &  zsoil,soiltyp,sigmaf,deltim,lat)
#endif
#endif
c
c  total runoff is composed of drainage into water table and
c  runoff at the surface and is accumulated in unit of meters
c
      do j = 1, LONS2                          
        runoff(j,lat) = runof(j) * dtf / 1000.         
     &                + runoff(j,lat)                            
        bgrun(j,lat) = drain(j) * dtf / 1000.        
     &                + bgrun(j,lat)
#ifdef ISOTOPE
        runoff1(j,lat) = rrsvi1(j,1)*runof(j) * dtf / 1000.         
     &                + runoff1(j,lat)                            
        bgrun1(j,lat) = rrsvi1(j,lsoil_)*drain(j) * dtf / 1000.        
     &                + bgrun1(j,lat)
        runoff2(j,lat) = rrsvi2(j,1)*runof(j) * dtf / 1000.         
     &                + runoff2(j,lat)                            
        bgrun2(j,lat) = rrsvi2(j,lsoil_)*drain(j) * dtf / 1000.        
     &                + bgrun2(j,lat)
#endif
        snowmelt(j,lat) = snowmelt(j,lat) + 
#ifdef HYDRO
     1                    snowmt(j)*dtf*hfus*rhoh2o
#else
#ifdef NOAHYDRO
     1                    snowmt(j)*hfus*rhoh2o
#else
     1                    snowmt(j)*deltim*hfus*rhoh2o
#endif
#endif
#ifdef RIVER
c  for river discharge calculation, snapshots of total runoff are stored
c    (kei)
        trunof(j,lat)=runof(j)+drain(j)
#ifdef ISOTOPE
        trunof1(j,lat)=rrsvi1(j,1)*runof(j)+
     $       rrsvi1(j,lsoil_)*drain(j)
        trunof2(j,lat)=rrsvi2(j,1)*runof(j)+
     $       rrsvi2(j,lsoil_)*drain(j)
#endif
#endif
      enddo                                          
c
c  return updated smsoil and stsoil to global arrays
c
      do k = 1, lsoil_
        do j = 1, LONS2
          smc(j,lat,k) = smsoil(j,k)
#ifndef VIC
          stc(j,lat,k) = stsoil(j,k)
#endif
#ifdef NOALSM1
          slc(j,lat,k) = slsoil(j,k)
#endif
#ifdef VICLSM1
          sic(j,lat,k) = sicld(j,k)
#endif
        enddo
      enddo
#ifdef VICLSM1
c
      do k = 1, nsoil_
         do j = 1, LONS2
            stc(j,lat,k) = stsoil(j,k)
         enddo
      enddo
      do j = 1, LONS2
         lstsn(j,lat) = lstsnld(j)
      enddo
#endif
c
      endif
c
c  zonal diagnostics
c
#ifdef DGZ
#ifndef SCM
      if(lastep) then
      secphy=shour
      secrad=max(shour,3600.*dtswav)
      call znlavb(lat,LONS2,LONF2S,
     1        levs_,secphy,secrad,
     &        rbs2(lat),sihyb,slhyb,delhyb,psexp,
     $        geshem(1,lat),bengsh(1,lat),
     &        dusfc(1,lat),dvsfc(1,lat),dtsfc(1,lat),dqsfc(1,lat),
     &        dugwd(1,lat),dvgwd(1,lat),
     &        slmsk(1,lat),sheleg(1,lat),tsea(1,lat),
     &        smc(1,lat,lsoil_),stc(1,lat,1),stc(1,lat,lsoil_),
     &        tg3(1,lat),zorl(1,lat),ep(1,lat),cldwrk(1,lat),
     &        dlwsfc(1,lat),ulwsfc(1,lat),
     &        xgrs,gt0,gq0,plamgr,pphigr,gu0,gv0)
      endif
#endif
#endif
c
c  grid point diagnostics
c
#ifdef KEN
      call kenput(lat,rcl,
     &        slmsk(1,lat),psexp,tg3(1,lat),sheleg(1,lat),radsl,dlwsf1,
#ifndef OSULSM1
     &        tsea(1,lat),qss,vegtype,gflx,zorl(1,lat),cd,cdq,
#else
     &        tsea(1,lat),qss,plantr(1,lat),gflx,zorl(1,lat),cd,cdq,
#endif
     &        rnet,hflx,stsoil,
     &        canopy(1,lat),drain,smsoil,runof,cld1d,
     &        u10m(1,lat),v10m(1,lat),t2m(1,lat),q2m(1,lat),
     &        hpbl(1,lat),gamt,gamq,
     &        dqsfc1,dtsfc1,dusfc1,dvsfc1,
     &        dusfcg,dvsfcg,
     &        rainc,rainl,
     &        gu0(1,1),gv0(1,1),
     &        gt0(1,1),gq0(1,1),ntotal_,
     &        hsw,hlw(1,1,lat),vvel,qcicnv,qrscnv,
     &        snowmt,snowev,snowfl)
#endif
c
      do j=1,LONS2
        pwat(j,lat)=0.
        qu00(j)=0.
        qv00(j)=0.
#ifdef ISOTOPE
        pwati1(j,lat)=0.
        pwati2(j,lat)=0.
#endif
      enddo
c
c  restore dry temp back to virtual temperature
c
#undef VTEMP
#ifdef RSM
#define VTEMP
#endif
#ifdef PSPLIT
#define VTEMP
#endif
#ifdef VTEMP
      do k=1,levs_
        do j=1,LONS2
          w2(j,k)=1.+fv*qgrs(j,k)
          work1(j) = 0.
          if(ncloud.gt.0) then
             do ic = icloud,ntotal_
               kc = (ic-1)*levs_ + k
#ifndef NISLQ
               work1(j) = max(qgrs(j,kc),qmin) + work1(j)
#else
               work1(j) = max(qgrs(j,kc),qmin) + work1(j)
#endif
             enddo
             w2(j,k)=w2(j,k) - work1(j)
          endif
        enddo
      enddo
      do k=1,levs_
        do j=1,LONS2
          tgrs(j,k)=tgrs(j,k)*w2(j,k)
        enddo
      enddo
#endif
      do k=1,levs_
        do j=1,LONS2
          w2(j,k)=1.+fv*gq0(j,k)
          work1(j) = 0.
          if(ncloud.gt.0) then
             do ic = icloud,ntotal_
               kc = (ic-1)*levs_ + k
#ifndef NISLQ
               work1(j) = max(gq0(j,kc),qmin) + work1(j)
#else
               work1(j) = max(gq0(j,kc),qmin) + work1(j)
#endif
             enddo
             w2(j,k)=w2(j,k) - work1(j)
          endif
        enddo
      enddo
      do k=1,levs_
        do j=1,LONS2
          gt0(j,k)=gt0(j,k)*w2(j,k)
        enddo
      enddo
c
      do k=1,levs_
        do j=1,LONS2
          pwat(j,lat)=pwat(j,lat)+delhyb(j,k)*gq0(j,k)
          qu00(j)=qu00(j)+delhyb(j,k)*gq0(j,k)*gu0(j,k)*rcs
          qv00(j)=qv00(j)+delhyb(j,k)*gq0(j,k)*gv0(j,k)*rcs
#ifdef ISOTOPE
          pwati1(j,lat)=pwati1(j,lat)+delhyb(j,k)*gq0(j,k+2*levs_)
          pwati2(j,lat)=pwati2(j,lat)+delhyb(j,k)*gq0(j,k+3*levs_)
#endif
        enddo
      enddo
      do j=1,LONS2
        pwat(j,lat)=pwat(j,lat)*psexp(j)*(1.e3/g_)
#ifdef HYDRO
        qull(j,lat)=qull(j,lat)+qu00(j)*psexp(j)*(1.e3/g_)*dtf
        qvll(j,lat)=qvll(j,lat)+qv00(j)*psexp(j)*(1.e3/g_)*dtf
#else
        qull(j,lat)=qull(j,lat)+qu00(j)*psexp(j)*(1.e3/g_)*deltim
        qvll(j,lat)=qvll(j,lat)+qv00(j)*psexp(j)*(1.e3/g_)*deltim
#endif
#ifdef ISOTOPE
        pwati1(j,lat)=pwati1(j,lat)*psexp(j)*(1.e3/g_)
        pwati2(j,lat)=pwati2(j,lat)*psexp(j)*(1.e3/g_)
        if ((pwati1(j,lat)/pwat(j,lat).ge.0.).and.
     $       (pwati1(j,lat)/pwat(j,lat).lt.1.04)) then
        else
           print*,'CAUTION! ISOPWAT IS BLOWING!!'
           print'(3i3,3e15.7)',mype,j,lat,
     $          pwati1(j,lat),pwat(j,lat),pwati1(j,lat)/pwat(j,lat)
           if (pwati1(j,lat).gt.-1.e10.and.pwati1(j,lat).lt.1.e10) then
              continue
           else
              print*,'Error! Abort'
              do k=1,levs_
                 print*,k,gq0(j,k+2*levs_)
              enddo
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
        endif
#endif
      enddo
#ifndef RSM
#ifndef SCM
c
c  quasi-pressure surface horizontal diffusion correction for t and q
c
#ifndef NISLQ
      call gridiff(LONS2,LONF2S,
     1             levs_,deltim,slhyb,pslap,gt0,gq0,ntotal_)
#endif

#ifdef NISLQ_PHYS
!
! replace negative with zero
!
      do k = 1,levh_
        do i = its,ite
          if(gq0(i,k).lt.0.) gq0(i,k) = 0.
        enddo
      enddo
#endif


#endif
#endif
#ifdef ISOTOPE
*** isotope correction 2***
      do k=1,levs_
        do j=1,LONS2
          gq0(j,k+levs_)=0.D0 !! dummy
          if (gq0(j,k).le.0.D0) then
            gq0(j,k)=0.D0
            gq0(j,k+2*levs_)=0.D0
            gq0(j,k+3*levs_)=0.D0
          else
            gq0(j,k+2*levs_)=max(gq0(j,k+2*levs_),0.D0)  ! tracer
            gq0(j,k+3*levs_)=max(gq0(j,k+3*levs_),0.D0)  ! tracer
          endif
c
#ifndef NOFRAC
          if ((gq0(j,k).le.1.D-5).and.(gq0(j,k).gt.0.D0).and.
     $           ((gq0(j,k+2*levs_)/gq0(j,k).gt.1.05D0).or.
     $            (gq0(j,k+3*levs_)/gq0(j,k).gt.1.40D0))) then
            alpha=1.-1./(-log10(gq0(j,k)/1.D-5)+1)**2.
            gq0(j,k+2*levs_)=gq0(j,k+2*levs_)+
     $           (gq0(j,k)-gq0(j,k+2*levs_))*alpha
            gq0(j,k+3*levs_)=gq0(j,k+3*levs_)+
     $           (gq0(j,k)-gq0(j,k+3*levs_))*alpha
          endif
#endif
        enddo
      enddo
c
      call isorsv(LONS2,LONF2S,dtf,
     &            rain,raini,
     $            smsoil,smsoil0,canopy(1,lat),sheleg(1,lat),
     $            drain,
     &            rrsvi1,rrsvi2)
c
      do k = 1, lsoil_
        do j = 1, LONS2
          rresv1(j,lat,k) = rrsvi1(j,k)
          rresv2(j,lat,k) = rrsvi2(j,k)
        enddo
      enddo
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(smsoil(1,1),LONS2,
     $            1,1,1,'smc1')
      call maxmin(smsoil(1,2),LONS2,
     $            1,1,1,'smc2')
      call maxmin(smsoil(1,3),LONS2,
     $            1,1,1,'smc3')
      call maxmin(smsoil(1,4),LONS2,
     $            1,1,1,'smc4')
      call maxmin(canopy(1,lat),LONS2,
     $            1,1,1,'canopy')
      call maxmin(sheleg(1,lat),LONS2,
     $            1,1,1,'sheleg')
      call maxmin(rain(1),LONS2,
     $            1,1,1,'rain')
      call maxmin(raini(1,1),LONS2,
     $            1,1,1,'rain1')
      call maxmin(raini(1,2),LONS2,
     $            1,1,1,'rain2')
      call maxmin(rrsvi1(1,1),LONS2,
     $            1,1,1,'rresv1')
      call maxmin(rrsvi2(1,1),LONS2,
     $            1,1,1,'rresv2')
#ifdef MP
      endif
#endif
#endif
#endif
c
#ifndef RSM
c
c  rcs2(lat) will be multiplied to ggu0 and ggv0 in gfidiu
c
      do k=1,levs_
        do j=1,LONS2
#ifndef PSPLIT
          ggu0(j,k)=gu0(j,k)
          ggv0(j,k)=gv0(j,k)
          ggt0(j,k)=gt0(j,k)
#else
          ggu0(j,k)=(gu0(j,k)-ugrs(j,k))/dt2
          ggv0(j,k)=(gv0(j,k)-vgrs(j,k))/dt2
          ggt0(j,k)=(gt0(j,k)-tgrs(j,k))/dt2
#endif
        enddo
      enddo
      do k=1,levh_
        do j=1,LONS2
#ifndef PSPLIT
          ggq0(j,k)=gq0(j,k)
#else
          ggq0(j,k)=(gq0(j,k)-qgrs(j,k))/dt2
#endif
        enddo
      enddo
#else
#ifdef PSPLIT
      do k=1,levs_
        do j=1,LONS2
c
c  rcs2(lat) will be multiplied to ggu0 and ggv0 in gfidiu
c
          ggu0(j,k)=(gu0(j,k)-ugrs(j,k))/dt2
          ggv0(j,k)=(gv0(j,k)-vgrs(j,k))/dt2
          ggt0(j,k)=(gt0(j,k)-tgrs(j,k))/dt2
        enddo
      enddo
      do k=1,levh_
        do j=1,LONS2
          ggq0(j,k)=(gq0(j,k)-qgrs(j,k))/dt2
        enddo
      enddo
#endif
#endif
#ifdef NFDR
c
c  store 2pm values 
c
      dt2hr=dt2/3600.
      do j=1,LONS2
#ifdef RSM
        rloctime=solhr+rlon(j,lat)/0.01745329/15.
#else
        rloctime=solhr+xlon(j,lat)/0.01745329/15.
#endif
        if(rloctime.lt. 0.) rloctime=rloctime+24.
        if(rloctime.ge.24.) rloctime=mod(rloctime,24.)
        if(rloctime-dt2hr.le.twopm.and.rloctime+dt2hr.gt.twopm) then
           t2m2pm(j,lat)=t2m(j,lat)
           rh2m2pm(j,lat)=rh2m(j)
           wind2pm(j,lat)=sqrt(u10m(j,lat)**2+v10m(j,lat)**2)
        endif
      enddo
#endif
c
      return
      end
