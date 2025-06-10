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
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define RD rd
#define RV rv
#define RERTH rerth
#define TWOJ1 twoj1
#define PI pi
      subroutine gloopn
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gloop  computes dynamic non-linear tendency terms
c of temp. div. ln(ps)
c computes predicted values of vorticity and moisture
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-13
c
c abstract:
c   program nudges the Fourier coefficients of forecasted run ("fcst data")
c   with those of reanalysis (or other external) run (called "base data").
c   program starts with the conversion of spectrum base data to 
c   fourier series. this conversion is taken place in every 6 hours
c   (nesthour). the base data is linearly interpolated for the time
c   in between. the spectra of fcst data are also converted to
c   fourier series.
c   two kinds of fourier series are nudged with a weighting constant as,
c    fnew=(1/(a+1))*(fold+a*fbase)
c   for those with longer wavelength than lnest (e.g. 1000km)
c   (for those with shorter, a=0).
c   a can be variable, too (attenuating function).
c   all of the variables, t,q,u,v,ps,gz can be nudged, but practically,
c   ps,u,v and t should be nudged. instead, zonal q (the first coefficient
c   of the fourier series) is replacable to that of the base data.
c
c default setting:
c   u,v,t,ps : nudging
c   q        : nothing (used to be zonal average replace)
c   n-scale  : 2000 km
c   n-weight : 0.9 (indaa=0)
c              used to be attenuating function (indaa=1)
c
c program history log:
c   07-07-24  kei yoshimura             modified from gloopa
c
c usage:    call gloopa
c
c attributes:
c   language: fortran 77
c
#include <paramodel.h>
#ifdef MP
#define LEVSS levsp_
#define LEVHS levhp_


c
c  begin twoloop(comfibm)
c
c  version with stacked transforms
c
c
c  begin comfspec
c
       common /comfspec/ idate(4),relvor(levs_),absvor(levs_),
     1 eps(lnuv_),epsi(lnuv_),gz(lnt22p_)
       common /comfspec/
     *    zem    (lnt22p_,levs_),
     *    dim    (lnt22p_,levs_),
     *    tem    (lnt22p_,levs_),
#ifndef NISLQ
     *    rm     (lnt22p_,levh_),
#endif
     *    qm     (lnt22p_)
       common /comfspec/
     *    ze     (lnt22p_,levs_),
     *    di     (lnt22p_,levs_),
     *    te     (lnt22p_,levs_),
#ifndef NISLQ
     *    rq     (lnt22p_,levh_),
#endif
     *    uln    (lnt22p_,levs_),
     *    vln    (lnt22p_,levs_),
     *    dpdphi (lnt22p_),
     *    dpdlam (lnt22p_),
     *    q      (lnt22p_),
     *    qlap   (lnt22p_)
       common /comfspec/
     *    z      (lnt22p_),
     *    uu     (lnt22p_,levs_),
     *    vv     (lnt22p_,levs_),
     *    y      (lnt22p_,levs_),
#ifndef NISLQ
     *    rt     (lnt22p_,levh_),
#endif
     *    x      (lnt22p_,levs_),
     *    w      (lnt22p_,levs_)
#ifdef NISLQ
       common /comfspec/ rm (lnt22p_,levh_)
       common /comfspec/ rq (lnt22p_,levh_)
       common /comfspec/ rt (lnt22p_,levh_)
#endif
c
c  begin comfgrid
c
       common /comfgrid/
     * colrad(latg2_),wgt(latg2_),wgtcs(latg2_),rcs2(latg2_),
     * sinlat(latg2_),
     * colrab(latg2p_),wgb(latg2p_),wgbcs(latg2p_),rbs2(latg2p_),
     * sinlab(lonf2p_,latg2p_),coslab(lonf2p_,latg2p_)
c
c
      logical lastep
c
c  begin comfver
c
      common/comfver/am(levs_,levs_),hm(levs_,levs_),tm(levs_,levs_),
     1               bm(levs_,levs_),cm(levs_,levs_),spdmax(levs_),
     2 si(levp1_),sl(levs_),del(levs_),rdel2(levs_),rmsdot(levm1_),
     3 ci(levp1_),cl(levs_),tov(levs_),gv(levs_),sv(levs_),rpi(levm1_),
     4 p1(levs_),p2(levs_), h1(levs_),   h2(levs_),rpirec(levm1_),
     5   thour,deltim,sl1,z00,fhour,shour,dtcvav,dtpost,
     6   filta,filtb,dk,tk,percut,dtswav,dtlwav,avprs0,
     7   cowave,dtwave
      common/comfveri/ kdt,jdt,inistp,limlow,maxstp,numsum,nummax
      common/comfverl/ stepone
      logical stepone
c
      common/comgpln/ qtt(lln2p_,latg2_),qvv(lcapp_,latg2_),
     1                qdd(lln2p_,latg2_),qww(lln2p_,latg2_)
c
      common/combitr/snnp1(lnt2_),runid,usrid
      common/combiti/ndex(lnt2_),ifin,icen,igen,icen2,
     1              ienst,iensi,ienss
      logical iope
      common/combitl/iope
      character*8 lab,labs
      common/combitc/lab(4),labs(4)
c
      common /radiag/ fluxr(lonf2p_,latg2p_,26)
      common /radiag/ cvavg(lonf2p_,latg2p_)
      parameter(ncpus=ncpus_)
c  sof include
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
      common/comdhc/spdlat(levs_,latg2p_)
c ...................................................
c local array
c................................................
c
      include 'mpif.h'
      integer npes,ncol,nrow,mype,master,msgtag
     *       ,myrow,mycol,comm_row,comm_column
     *       ,levstr,levlen,lonstr,lonlen,latstr,latlen
     *       ,lwvstr,lwvlen,lntstr,lntlen,lnpstr,lnplen
     *       ,lwvdef,latdef
      common /commpi/ npes,ncol,nrow,mype,master,msgtag
     *               ,myrow,mycol,comm_row,comm_column
     *               ,levstr(0:npes_-1),levlen(0:npes_-1)
     *               ,lonstr(0:npes_-1),lonlen(0:npes_-1)
     *               ,latstr(0:npes_-1),latlen(0:npes_-1)
     *               ,lwvstr(0:npes_-1),lwvlen(0:npes_-1)
     *               ,lntstr(0:npes_-1),lntlen(0:npes_-1)
     *               ,lnpstr(0:npes_-1),lnplen(0:npes_-1)
     *               ,lwvdef(jcap_+1)  ,latdef(latg_/2)
c
c   begin compspec
c
        common /compspec/
     *    zea    (lln22p_,levsp_),
     *    dia    (lln22p_,levsp_),
     *    tea    (lln22p_,levsp_),
#ifndef NISLQ
     *    rqa    (lln22p_,levhp_),
#endif
     *    ulna   (lln22p_,levsp_),
     *    vlna   (lln22p_,levsp_),
     *    dpdphia(lln22p_),
     *    dpdlama(lln22p_),
     *    qa     (lln22p_),
     *    qlapa  (lln22p_)
       common /compspec/
     *    za     (lln22p_),
     *    uua    (lln22p_,levsp_),
     *    vva    (lln22p_,levsp_),
     *    ya     (lln22p_,levsp_),
#ifndef NISLQ
     *    rta    (lln22p_,levhp_),
#endif
     *    xa     (lln22p_,levsp_),
     *    wa     (lln22p_,levsp_)
#ifdef NISLQ
        common /compspec/ rqa (lln22p_,levhp_)
        common /compspec/ rta (lln22p_,levhp_)
#endif
c
      dimension
     1 syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_),
     3 anf(lonf22_,lotas,latg2p_)
c
      integer lcapd,lonfd,lonfdp,lonfds,lcapdp
      common /comreduce/
     *                lcapd (latg_/2) ! regular grid in lat
     *               ,lonfd (latg_/2) ! regular grid in lat
     *               ,lcapdp(latg_/2,0:npes_-1) !regular grid in lat
     *               ,lonfdp(latg2p_,0:npes_-1) !symmetry spread grid in lat
     *               ,lonfds(latg2p_,0:npes_-1) !symmetry spread grid in lat
c ........................................................
      dimension scs2(jcap1p_),
     2 syn(lcap22p_,lotss,latg2_),syntop(2,jcap1p_,lotsts),
     4 anl(lcap22p_,lotas,latg2_),anltop(2,jcap1p_,lotats),
     5 flp(2,jcap1p_,lotas,latg2_),flm(2,jcap1p_,lotas,latg2_)
cccKEIccc
      character*20 bfname,bfname2
      real fnxx,sshour,nestsec,sshour2,fnxx2,fnxb
      data sshour,sshour2 /0.,0./
      data fnxb/-1./
      real fhour1,fhour2,nestsec2
      integer nchar,idate1(4),idate2(4),nbase1,nbase2
      data nbase1,nbase2/91,92/
      data nchar/15/
      integer ifirst
      data ifirst/0/
      integer m,stm,edm,ii,jj1,jj2,jlim(latg2p_)
      real nlat(latg2_),nlatb(latg2p_)
      save jlim,nlat,nlatb,ifirst,sshour,fnxb,bfname2
      real lnest,nesthour, aa   !! nudging scale, interval, weight
      data lnest,nesthour, aa
     $    /2000.e3,  6.  ,0.9/
      real aa2
      integer jnest,indaa       !! basedata resolution, use of AF for aa
      data    jnest,indaa       !! indaa=0: use constant value.
     $       /   62,    0/
      real
     *    gz1     (lnt22p_),
     *     q1     (lnt22p_),
     *    ze1     (lnt22p_,levs_),
     *    di1     (lnt22p_,levs_),
     *    te1     (lnt22p_,levs_),
     *    rq1     (lnt22p_,levh_),
     *    z001
      real
     *    gz2     (lnt22p_),
     *     q2     (lnt22p_),
     *    ze2     (lnt22p_,levs_),
     *    di2     (lnt22p_,levs_),
     *    te2     (lnt22p_,levs_),
     *    rq2     (lnt22p_,levh_),
     *    z002
      real
     *    bzea1    (lln22p_,levsp_),
     *    bdia1    (lln22p_,levsp_),
     *    btea1    (lln22p_,levsp_),
     *    brqa1    (lln22p_,levhp_),
     *    bulna1   (lln22p_,levsp_),
     *    bvlna1   (lln22p_,levsp_),
     *    bdpdphia1(lln22p_),
     *    bdpdlama1(lln22p_),
     *    bqa1     (lln22p_),
     *    bqlapa1  (lln22p_)
      real
     *    bzea2    (lln22p_,levsp_),
     *    bdia2    (lln22p_,levsp_),
     *    btea2    (lln22p_,levsp_),
     *    brqa2    (lln22p_,levhp_),
     *    bulna2   (lln22p_,levsp_),
     *    bvlna2   (lln22p_,levsp_),
     *    bdpdphia2(lln22p_),
     *    bdpdlama2(lln22p_),
     *    bqa2     (lln22p_),
     *    bqlapa2  (lln22p_)
c
      dimension byf(lonf22_,lotss,latg2p_)
      dimension byf1(lonf22_ ,lotss,latg2p_),
     2 byn1(lcap22p_,lotss,latg2_),byntop1(2,jcap1p_,lotsts)
      dimension byf2(lonf22_ ,lotss,latg2p_),
     2 byn2(lcap22p_,lotss,latg2_),byntop2(2,jcap1p_,lotsts)
      save byf1,byf2
      real ditmp(lnt22p_,levs_)
c
      dimension ek(lonf22p_,levs_)      
c
c#undef latg2_
c
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
      jstr=latstr(mype)
      jend=latstr(mype)+latlen(mype)-1
      lons2=lonlen(mype)*2
      lats2=latlen(mype)
c
c=========================================================================
c GSM : calculate dynamic tendency terms
c=========================================================================
c
cccKEIccc
c  nudging by base field
c   read base
c
      nestsec=nesthour*60.*60.
      sshour2=sshour
      if (ifirst.eq.0) then
        if (iope) then
           print*,'gloopn, lnest,aa,indaa',lnest,aa,indaa
        endif
        if (shour.gt.1800.) then
          sshour=shour+(thour-1.)*60.*60.
          sshour2=(thour-1.)*60.*60.
        else
          sshour=shour+(thour)*60.*60.
          sshour2=(thour)*60.*60.
        endif
        do j=1,latg2_
** only tropics (latitudinal weighting)***
          if (colrad(j).lt.0.866) then
            nlat(j)=1.
          else
            nlat(j)=0.
         endif
        enddo
c
        do j=1,latlen(mype)
          jj=latdef(latstr(mype)+j-1)
          nlatb(j)=nlat(jj)
** with latitudinal weighting**
c          jlim(j)=min(int((2*pi_*rerth_*sin(colrad(jj))
c     $          /lnest)*nlatb(j)),jnest)
          jlim(j)=min(int(2*pi_*rerth_*sin(colrad(jj))
     $          /lnest),jnest)
        enddo
        lat1=jstr
        lat2=jend
        latdon=jstr-1
        do lat=lat1,lat2
          lan=lat-latdon
          if (iope) then
          print 104,'JJ',mype,lat,jlim(lan),
     $         nlatb(lan),colrab(lan)
          endif
        enddo
c
        fnxx=int(sshour2/nestsec)*nesthour
        if (iope) then
           call basefnam(idate(4),idate(2),idate(3),idate(1),
     1          fnxx,
     2          bfname,nchar)
           close(nbase1)
           open(unit=nbase1,file=bfname(1:nchar),status='unknown',
     1          form='unformatted')
        endif
c
c        call rdsig2(nbase1,fhour1,idate1,gz1,q1,te1,di1,ze1,rq1,
        call rdsig(nbase1,fhour1,idate1,gz1,q1,te1,di1,ze1,rq1,
     $       sl,si,z001)
*** without divergence component ***
c        do j=1,lnt22p_
c          do k=1,levs_
c            di1(j,k)=0.
c          enddo
c        enddo
c
      else
        sshour=sshour+deltim
      endif
 104  format (a5,3i4,f6.2,f15.7)
      fnxx=int(sshour2/nestsec)*nesthour
      fnxx2=fnxx+nesthour
c   read basefiles (sigma-level spectral coef.)
      if (fnxb.ne.fnxx) then
c
        if (iope) then
c
        call basefnam(idate(4),idate(2),idate(3),idate(1),
     1       fnxx2,
     2       bfname2,nchar)
 113    format(a5,2f8.2,2f10.0,' ',a15,2f8.2)
        close(nbase2)
        open(unit=nbase2,file=bfname2(1:nchar),status='unknown',
     1       form='unformatted')
c
        endif
c
c        call rdsig2(nbase2,fhour2,idate2,gz2,q2,te2,di2,ze2,rq2,
        call rdsig(nbase2,fhour2,idate2,gz2,q2,te2,di2,ze2,rq2,
     $        sl,si,z002)
*** without divergence component ***
c        do j=1,lnt22p_
c          do k=1,levs_
c            di2(j,k)=0.
c          enddo
c        enddo
      endif
      if (iope) then
         print 113,'Kei',fnxx,fnxx2,sshour,sshour2,bfname2,thour,fnxb
      endif
c
*** without divergence component ***
c      do j=1,lnt22p_
c         do k=1,levs_
c            ditmp(j,k)=di(j,k)
c            di(j,k)=0.
c         enddo
c      enddo
c
#ifndef NISLQ
      call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,3+ntotal_)
#else
      call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,3)
#endif
      call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
      call dellnp(qa,dpdphia,syntop(1,1,2*levsp_+1),dpdlama,
     1            llstr,llens,lwvdef)
      call deldfsp(qa,qlapa,llstr,llens,lwvdef)
      call dzuvle(dia,zea,ulna,vlna,syntop(1,1,1),
     1            syntop(1,1,levsp_+1),llstr,llens,lwvdef)
c
cc for base cc
      if (fnxb.ne.fnxx) then
c
      if (ifirst.eq.0) then
c
      call mpnn2nk(ze1 ,lnt22p_,levs_, bzea1,lln22p_,levsp_,1)
      call mpnn2nk(di1 ,lnt22p_,levs_, bdia1,lln22p_,levsp_,1)
      call mpnn2nk(te1 ,lnt22p_,levs_, btea1,lln22p_,levsp_,1)
      call mpnn2nk(rq1 ,lnt22p_,levs_, brqa1,lln22p_,levhp_,1)
      call mpnn2n(q1 ,lnt22p_, bqa1,lln22p_,1)
      call dellnp(bqa1,bdpdphia1,byntop1(1,1,2*levsp_+1),bdpdlama1,
     1            llstr,llens,lwvdef)
      call deldfsp(bqa1,bqlapa1,llstr,llens,lwvdef)
      call dzuvle(bdia1,bzea1,bulna1,bvlna1,byntop1(1,1,1),
     1            byntop1(1,1,levsp_+1),llstr,llens,lwvdef)
      endif
c
      call mpnn2nk(ze2 ,lnt22p_,levs_, bzea2,lln22p_,levsp_,1)
      call mpnn2nk(di2 ,lnt22p_,levs_, bdia2,lln22p_,levsp_,1)
      call mpnn2nk(te2 ,lnt22p_,levs_, btea2,lln22p_,levsp_,1)
      call mpnn2nk(rq2 ,lnt22p_,levs_, brqa2,lln22p_,levhp_,1)
      call mpnn2n(q2 ,lnt22p_, bqa2,lln22p_,1)
      call dellnp(bqa2,bdpdphia2,byntop2(1,1,2*levsp_+1),bdpdlama2,
     1            llstr,llens,lwvdef)
      call deldfsp(bqa2,bqlapa2,llstr,llens,lwvdef)
      call dzuvle(bdia2,bzea2,bulna2,bvlna2,byntop2(1,1,1),
     1            byntop2(1,1,levsp_+1),llstr,llens,lwvdef)
c
      endif
c
      do n=1,2
        do j=1,jcap1p_
          do l=1,lotats
            anltop(n,j,l)=0.0
          enddo
        enddo
      enddo
c
      lat1=1
      lat2=latg2_
      latdon=0
c first lat loop
!$omp parallel do private(lat,lan,llensd)
c
      do 1000 lat=lat1,lat2
        lan=lat-latdon
        llensd=lcapdp(lat,mype)
        call sums2i(zea,syn(1,1,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,lotss)
        call sumtop(syn(1,kstbs,lan),syntop,qvv(1,lat),
     &              llstr,llensd,lwvdef,lotsts)
cc for base cc
        if (fnxb.ne.fnxx) then
c
        if (ifirst.eq.0) then
c
        call sums2i(bzea1,byn1(1,kszs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bdia1,byn1(1,ksds,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(btea1,byn1(1,ksts,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
#ifndef NISLQ
        call sums2i(brqa1,byn1(1,ksrs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levhp_)
#endif
        call sums2i(bulna1,byn1(1,ksus,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bvlna1,byn1(1,ksvs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bqa1,byn1(1,ksps,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,1)
        call sumtop(byn1(1,kstbs,lan),byntop1,qvv(1,lat),
     &              llstr,llensd,lwvdef,lotsts)
        endif
c
        call sums2i(bzea2,byn2(1,kszs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bdia2,byn2(1,ksds,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(btea2,byn2(1,ksts,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
#ifndef NISLQ
        call sums2i(brqa2,byn2(1,ksrs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levhp_)
#endif
        call sums2i(bulna2,byn2(1,ksus,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bvlna2,byn2(1,ksvs,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,levsp_)
        call sums2i(bqa2,byn2(1,ksps,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,1)
        call sumtop(byn2(1,kstbs,lan),byntop2,qvv(1,lat),
     &              llstr,llensd,lwvdef,lotsts)
        endif
c       
1000  continue
c
c compute merid. derivs. of temp. and moisture using qdd.
c
      call mpnl2ny(syn,lcap22p_,latg2_,
     1             syf,lonf22_,latg2p_,lotss,1,lotss)
cc for base cc
      if (fnxb.ne.fnxx) then
c
      if (ifirst.eq.0) then
        ifirst=1
c
        call mpnl2ny(byn1,lcap22p_,latg2_,
     1             byf1,lonf22_,latg2p_,lotss,1,lotss)
c
      else
        do i=1,lonf22_
          do k=1,lotss
            do j=1,latg2p_
              byf1(i,k,j)=byf2(i,k,j)
            enddo
          enddo
        enddo
c
      endif
c
      call mpnl2ny(byn2,lcap22p_,latg2_,
     1             byf2,lonf22_,latg2p_,lotss,1,lotss)
      endif
c
      do i=1,lonf22_
        do k=1,lotss
          do j=1,latg2p_
            byf(i,k,j)
     $            =byf1(i,k,j)*(1.-(sshour-fnxx*3600.)/nestsec)
     $            +byf2(i,k,j)*    (sshour-fnxx*3600.)/nestsec
          enddo
        enddo
      enddo
c
cc nudging cc
      lat1=jstr
      lat2=jend
      latdon=jstr-1
c
      do lat=lat1,lat2
        lan=lat-latdon
        do n=1,2
          do i=1,(jlim(lan)+1)*2
            ii=(n-1)*(lonf22_-2)/2+i
            if (indaa.eq.0) then
              aa2=aa/(1.+aa)
            else
              aa2=(1.-real(int((i-1)/2)/max(real(jlim(lan)),1.)))**0.5
                                !! alpha infinity->0 sqrt (YK07)
c              aa2=(cos(int((i-1)/2)/max(jlim(lan),1)*pi_)+1.)/2.
                                !! alpha Infinity->0 cosine
            endif
                                !! aa2 is transformed in aa/(1+aa)
            do k=1,levsp_
              syf(ii,ksus+k-1,lan)
     $             =syf(ii,ksus+k-1,lan)*(1.-aa2)
     $             +byf(ii,ksus+k-1,lan)*(   aa2) !! u
              syf(ii,ksvs+k-1,lan)
     $             =syf(ii,ksvs+k-1,lan)*(1.-aa2)
     $             +byf(ii,ksvs+k-1,lan)*(   aa2) !! v
              syf(ii,ksts+k-1,lan)
     $             =syf(ii,ksts+k-1,lan)*(1.-aa2)
     $             +byf(ii,ksts+k-1,lan)*(   aa2) !! t
            enddo
            syf(ii,ksps,lan)
     $           =syf(ii,ksps,lan)*(1.-aa2)
     $           +byf(ii,ksps,lan)*(   aa2) !! q (sfP)
          enddo
c          do i=1,(jlim(lan)+1)*2 !! for same as UV-SSBC
          do i=1,2              !! for zonal mean nudging
            ii=(n-1)*(lonf22_-2)/2+i
c            do k=1,levsp_
c              syf(ii,ksts+k-1,lan)=byf(ii,ksts+k-1,lan) !! te(tmp)
c            enddo
#ifndef NISLQ
c            do k=1,levhp_
c              syf(ii,ksrs+k-1,lan)=byf(ii,ksrs+k-1,lan) !! rq(vpr)
c            enddo
#endif
c            syf(ii,ksps,lan)=byf(ii,ksps,lan)           !! q (sfP)
          enddo
        enddo
      enddo
c      if (mype.eq.master) then
c         do i=1,lonf22_
c            print 102,'syf2',i,ksps,syf(i,ksps,10),byf(i,ksps,10)
c         enddo
c         do i=1,lonf22_
c            print 102,'syf2',i,ksus,syf(i,ksus,10),byf(i,ksus,10)
c         enddo
c         do i=1,lonf22_
c            print 102,'syf2',i,ksvs,syf(i,ksvs,10),byf(i,ksvs,10)
c         enddo
c         do i=1,lonf22_
c            print 102,'syf2',i,ksts,syf(i,ksts,10),byf(i,ksts,10)
c         enddo
c         do i=1,lonf22_
c            print 102,'syf2',i,ksrs,syf(i,ksrs,10),byf(i,ksrs,10)
c         enddo
c      endif
 102  format (a5,2i5,2e15.7)
c
!$omp parallel do private(j,k)
      do k=1,levsp_
       do j=1,lln22p_
        uua(j,k)=0.0
        vva(j,k)=0.0
         ya(j,k)=0.0
       enddo
      enddo
#ifndef NISLQ
      do k=1,levhp_
       do j=1,lln22p_
        rta(j,k)=0.0
       enddo
      enddo
#endif
      do j=1,lln22p_
       za(j)=0. e 0
      enddo
c
!$omp parallel do private(i,j,k)
      do i=1,lonf22_
         do j=1,latg2p_
            anf(i,kaps,j)=syf(i,ksps,j)            !! q (sfP)
            do k=1,levsp_
               anf(i,kaus-1+k,j)=syf(i,ksus-1+k,j) !! u
               anf(i,kavs-1+k,j)=syf(i,ksvs-1+k,j) !! v
               anf(i,kats-1+k,j)=syf(i,ksts-1+k,j) !! te (tmp)
            enddo
#ifndef NISLQ
            do k=1,levhp_
               anf(i,kars-1+k,j)=syf(i,ksrs-1+k,j) !! rt (vpr)
            enddo
#endif
         enddo
      enddo
c
      call mpny2nl(anf,lonf22_ ,latg2p_,
     1             anl,lcap22p_,latg2_ ,lotas,kaps,lotas)
c
      lat1=1
      lat2=latg2_
      latdon=0
c
!$omp parallel do private(lat,lan,llensd)
c
      do 2500 lat=lat1,lat2
        lan=lat-latdon
        llensd=lcapdp(lat,mype)
c
        do k=1,levsp_*2
          kaa=kaus-1+k
          do j=1,llensd*2
            anl(       j,kaa,lan)=anl(       j,kaa,lan)*rcs2(lat)
            anl(lcapp_+j,kaa,lan)=anl(lcapp_+j,kaa,lan)*rcs2(lat)
          enddo
        enddo
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
        llensd=lcapdp(lat,mype)
c
        call fl2ip(flp(1,1,1,lan),flm(1,1,1,lan),za,qww(1,lat),
     1             llstr,llensd,lwvdef,lotas)
c
        call uvsums(flp(1,1,kaus,lan),flm(1,1,kaus,lan),
     1              flp(1,1,kavs,lan),flm(1,1,kavs,lan),
     2              anltop(1,1,1),anltop(1,1,levsp_+1),
     3              qvv(1,lat),wgt(lat),
     4              llstr,llensd,lwvdef,levsp_)
c
3000  continue
c
#ifdef DBG
ccc Check UV (Unnecessary) ccc
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
c
       do 1300 lat=lat1,lat2
         lcapf=lcapd(latdef(lat))
         lonff=lonfd(latdef(lat))
         lan=lat-latdon
         call ftilonf(syf(1,1,lan),syf(1,1,lan),lotss*2,
     1                lcapf,lonff,latdef(lat),1)
 1300  continue
c
#ifndef NISLQ
       call mpnk2nx(syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,
     &              5+ntotal_)
#else
       call mpnk2nx(syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,
     &              5)
#endif
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
c
      do 2000 lat=lat1,lat2
        lan=lat-latdon
        lonsd2=lonfdp(lan,mype)*2
c
        do k=1,levs_
          spdlat(k,lan)=0.
        enddo
c
      do 140 k=1,levs_
      do 140 j=1,lonsd2
      ek(j,k)=(grs(j,ksu-1+k,lan)*grs(j,ksu-1+k,lan)
     $          +grs(j,ksv-1+k,lan)*grs(j,ksv-1+k,lan))*rbs2(lan)
      if (ek(j,k) .gt. spdlat(k,lan))  spdlat(k,lan)=ek(j,k)
 140  continue
c
 2000 continue
c
      do k=1,levs_
        spdmax(k) = 0.0
        do lat=1,lats2
          spdmax(k)=max(spdmax(k),spdlat(k,lat))
        enddo
        spdmax(k)=sqrt(spdmax(k))
      enddo
      call mpgetspd(spdmax)
      if( mype.eq.master ) then
      print 100,(spdmax(k),k=1,levs_)
100   format(' global checked speed maxima for all layers 2',
     &     :/' spdmx(01:10)=',10f5.0,:/' spdmx(11:20)=',10f5.0,
     &     :/' spdmx(21:30)=',10f5.0,:/' spdmx(31:40)=',10f5.0,
     &     :/' spdmx(41:50)=',10f5.0,:/' spdmx(51:60)=',10f5.0,
     &     :/' spdmx(61:70)=',10f5.0,:/' spdmx(71:80)=',10f5.0,
     &     :/' spdmx(81:90)=',10f5.0,:/' spdmx(91:00)=',10f5.0)
      endif
c
      if( mype.eq.master ) then
      do k=1,levs_
        if(spdmax(k).eq.0.) then
          print *,'run failure.  spdmax=0.'
          call abort
        endif
      enddo
      endif
c
ccc End of UV check ccc     
#endif
c
      call uvtodz(uua,vva,xa,wa,anltop(1,1,1),
     1            anltop(1,1,levsp_+1),llstr,llens,lwvdef)
c
      call mpnk2nn(ya,lln22p_,levsp_, te,lnt22p_,levs_,1)
#ifndef NISLQ
      call mpnk2nn(rta,lln22p_,levhp_, rq,lnt22p_,levh_,1)
#endif
      call mpnk2nn(xa,lln22p_,levsp_, di,lnt22p_,levs_,1)
*** without divergence component ***
c      do j=1,lnt22p_
c        do k=1,levs_
c          di(j,k)=ditmp(j,K)
c        enddo
c      enddo
      call mpnk2nn(wa,lln22p_,levsp_, ze,lnt22p_,levs_,1)
      call mpn2nn(za,lln22p_,q,lnt22p_,1)
c
      fnxb=fnxx
c
#endif
      return
      end
