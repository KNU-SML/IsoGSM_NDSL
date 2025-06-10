#define PI pi
#define LEVP1 levp1
#define BGF bgf
#define BORDER border
#define G g
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define IO io
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JO jo
#define KO ko
#define KT kt
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVS levs
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNTP lntp
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
#define NROW nrow
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define QMIN qmin
#define RD rd
#define RV rv
#include <define.h>
#ifndef RSM
#undef STDAMP_P
#endif
#ifndef HYBRID
      subroutine subpgb(fhour,idate,sl,si,
#else
      subroutine subpgb(fhour,idate,ak5,bk5,
#endif
     &                  fgz,fq,fte,fuu,fvv,frq,flat,flon,n1
#ifdef RSM
     &                 ,lbaseout,fm2,lskipgb)
#else
     &                  ,lskipgb)
#endif
#ifdef MRG_POST
c$$$  subprogram documentation block
c
c subprogram: subpgb        transforms a sigma to pressure grib
c   prgmmr: kanamitsu          date: 04-06-22
c
c abstract: transforms a sigma  spectral) file to pressure grib1
c
c usage:    call subpgb(fhour,idate,sl,si,
c    &                  fgz,fq,fte,fuu,fvv,frq,flat,flon,n1)
c   input arguments:
c     fhour        real forecast hour
c     idate        integer (4) date
c     si           real (levr+1) sigma interface values
c     sl           real (levr) sigma full level values
c     fgz
c     fq
c     fte
c     fuu
c     fvv
c     frq
c     flat
c     flon
c     n1
c
c subprograms called:
c   idsdef       set ids
c   sunprm       set parameters for sundry fields
c   rcof2grd     convert coef to grid
c   getrh        compute relative humidity
c   omega        compute vertical velocity
c   hydro        compute geopotential heights
c   sig2p        interpolate sigma to pressure
c   sundry       compute sundry fields
c   ptranw       quarterpack and transpose data
c   ptranr       unpack quarterpacked transposed data
c   gribit       create grib message
c   wryte        write data by bytes
c
c attributes:
c   language: fortran 90
c
c  define MPIGRIB will gather all the arrays to master processor then
c  convert to grib and writes on master processor.  this is most inefficient.
c
c  define PGBGATHER will conver to grib on each processor, then gather
c  gribbed arrays to master processor and then master processor writes.
c  this is second efficient, but i/o is serial.
c
c  undef MPIGRIB, undef PGBGATHER will do the coversion to grib on
c  each processor and then each processor writes.  this is most efficient
c  since computation and i/o are both parallel.
c
c  when MPIGRIB is defined, PGBGATHER should not be defined.
c  when PGBGATHER is defined, PGBGRIB should not be defined.
c
c$$$
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#define MPSYNALL rmpsynall
#else
#include <commpi.h>
#define MPSYNALL mpsynall
#endif
#endif
#ifndef RSM
#undef STDAMP_P
#endif
      parameter (l2=2)
      parameter (g= g_ ,rd= rd_ ,rv= rv_,fvirt=rv/rd-1.)
      parameter (gor2=g/rd*l2,eps=rd/rv,epsm1=rd/rv-1.)
      parameter (qmin=qmin_)
      parameter (nupa=11,nupt=4,nsun=29,nupas=8)
      parameter (lenpds=28,lengds=32)
c
      dimension idate(4),si(levs_+1),sl(levs_),po(ko_)
      dimension ids(255)
      dimension ipo(ko_),npo(ko_),pokpa(ko_)
      dimension npt1(kt_),npt2(kt_),iz(ko_)

#ifdef RSM
      logical lbaseout
#endif
c
#include <postplevs.h>
c
      dimension pot(255)
      data pot/255*0./
c
      parameter(rad=pi_/180.)
c
#ifdef RSM
#include <rscomloc.h>
#include <rscomltb.h>
#define IDIMF  igrd1_
#define JDIMF  jgrd1_
#define IDIMF2 igrd12_
#define JDIMF2 jgrd12_
#ifdef MP
#define LNWAVS lnwavp_
#define IDIM2 igrd12p_
#define JDIM2 jgrd12p_
#else
#define LNWAVS lnwav_
#define IDIM2 igrd1_
#define JDIM2 jgrd1_
#endif
#define IDIMO igrd1_
#define JDIMO jgrd1_

#else
#ifdef REDUCE_GRID
#include <comfgrid.h>
#include <comreduce.h>
#endif
#define IDIMF  lonf_
#define JDIMF  latg_
#define IDIMF2 lonf2_
#define JDIMF2 latg2_
#ifdef MP
#define LNWAVS lnt22p_
#define IDIM2 lonf2p_
#define JDIM2 latg2p_
#else
#define LNWAVS lnt22_
#define IDIM2 lonf2_
#define JDIM2 latg2_
#endif
#define IDIMO io_
#define JDIMO jo_
#endif
c
#ifdef COUPLE_RSM_AND_ROMS
#ifdef MP
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#define LNWAVS lnwavp_
#define LNGRDS lngrdp_
#else
#define IGRD12S igrd12_
#define JGRD12S jgrd12_
#define LNWAVS lnwav_
#define LNGRDS lngrd_
#endif
#ifndef RSM
#  ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#  else
#define LONF2S lonf2_
#define LATG2S latg2_
#  endif
#else
#  ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#  else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#  endif
#endif
#ifdef RSM
#include <comrflx.h>
#else
#include <comfphys.h>
#endif
#endif
c
      dimension fgz(LNWAVS),fq(LNWAVS),
     $          fte(LNWAVS,levs_),
     $          fuu(LNWAVS,levs_),fvv(LNWAVS,levs_),
     $          frq(LNWAVS,levh_)
#ifdef RSM
      dimension flat(IDIM2,JDIM2),flon(IDIM2,JDIM2)
#else
      dimension flat(1,JDIMF)
#endif
c
#ifdef RSM
#ifdef OSULSM1
#ifdef ISOTOPE
      parameter(nflx=83+2*lsoil_)
#else
      parameter(nflx=59+2*lsoil_)
#endif
#endif
#ifdef OSULSM2
#ifdef ISOTOPE
      parameter(nflx=83+2*lsoil_)
#else
      parameter(nflx=59+2*lsoil_)
#endif
#endif
#ifdef NOALSM1
#ifdef NOAHYDRO
#ifdef ISOTOPE
      parameter(nflx=84+2*lsoil_)
#else
      parameter(nflx=60+2*lsoil_)
#endif
#else
#ifdef ISOTOPE
      parameter(nflx=83+2*lsoil_)
#else
      parameter(nflx=59+2*lsoil_)
#endif
#endif
#endif
#ifdef VICLSM1
#ifdef ISOTOPE
      parameter(nflx=83+lsoil_+nsoil_)
#else
      parameter(nflx=59+lsoil_+nsoil_)
#endif
#endif
#else
      parameter(nflx=0)
#endif
c
      parameter(npgblen=nupa*ko_+nupt*kt_+nsun)
      parameter(ntlen=npgblen+nflx)
      dimension oxs(IDIM2,levs_),osxs(IDIM2,levs_)
      dimension rxs(IDIM2,levs_),qsxs(IDIM2,levs_)
      dimension zxs(IDIM2,levs_),zxi (IDIM2,levs_)
      dimension qxp(IDIM2,kt_)
      dimension fxp(IDIM2,ntlen)
c
      dimension ipu(ntlen)
      dimension itl(ntlen)
      dimension il1(ntlen)
      dimension il2(ntlen)
      dimension mxbit(ntlen)
      dimension ifhr(ntlen)
      dimension ithr(ntlen)
      dimension iavg(ntlen)
c
      parameter(nuptlen=nupt*kt_)
      logical lskipgb
      logical nsel(ntlen),lsel
      data  nsel
c              hgtprs     ugrdprs    vgrdprs    tmpprs     vvel
     $       / ko_*.true.,ko_*.true.,ko_*.true.,ko_*.true.,ko_*.true.,
c              rhprs      spfhprs    absv
#ifdef ES
     $          ko_*.true.,ko_*.true.,ko_*.false.,
#else
     $          ko_*.true.,ko_*.true.,ko_*.true.,
#endif
c              ozone
#if ( _ntrac_ == 2 )
     $          ko_*.true.,
#else
     $          ko_*.false.,
#endif
c              cldwat 
#if ( _ncldg_ >= 1 )
     $          ko_*.true.,
#else
#ifndef CO2PROG
#ifdef ISOTOPE
     $          ko_*.true.,
#else
     $          ko_*.false.,
#endif
#else
     $          ko_*.true.,
#endif
#endif
c              prcpwat
#if ( _ncldg_ >= 2 )
     $          ko_*.true.,
#else
#ifdef ISOTOPE
     $          ko_*.true.,
#else
     $          ko_*.false.,
#endif
#endif
c              dummy      
     $          nuptlen*.false.,
c              pressfc  ptend    pwat     rhclm     tmptrp   prestrp
#ifdef ES
     $         .true.,  .true.,  .true.,  .false.,  .false., .false.,
#else
     $         .true.,  .true.,  .true.,   .true.,  .true.,  .true.,
#endif
c              ugrdtrp  vgrdtrp  vssh     lftx     lftxb
#ifdef ES
     $         .false., .false., .false., .true.,  .true.,
#else
     $         .true.,  .true.,  .true.,  .true.,  .true.,
#endif
c              tmpmwl   presmwl  ugrdmwl  vgrdmwl   hgtsfc   prmsl
#ifdef ES
     $         .false., .false.,  .false., .false., .true.,  .true.,
#else
     $         .true.,  .true.,   .true.,  .true.,  .true.,  .true.,
#endif
c              dummy  
     $           10*.false.,
c              cldwat
#if ( _ncldg_ >= 1 )
     $           .true.,
#else
#ifndef CO2PROG
#ifdef ISOTOPE
     $           .true.,
#else
     $           .false.,
#endif
#else
     $           .true.,
#endif
#endif
c              ttlo3 
#if ( _ncldg_ >= 2 )
     $           .true. 
#else 
#ifdef ISOTOPE
     $           .true. 
#else
     $           .false.
#endif
#endif
#ifdef RSM
c              uflx     vflx     shtfl    lhtfl     tmpsfc
     $        ,.true.,  .true.,  .true.,  .true.,   .true.,
c              soilw             tmpdlr
#ifndef VIC
     $         lsoil_*.true.,    lsoil_*.true.,
#else
     $         lsoil_*.true.,    nsoil_*.true.,
#endif
c              weasd    dlwrfsfc ulwrfsfc ulwrftoa
     $         .true.,  .true.,  .true.,  .true.,
c              uswrftoa uswrfsfc dswrfsfc tcdchcl   preshct  preshcb
#ifdef ES
     $         .true.,  .true.,  .true.,  .true.,   .false., .false.,
#else
     $         .true.,  .true.,  .true.,  .true.,   .true.,  .true.,
#endif
c              tmphct   tcdcmcl  presmct  presmcb   tmpmct   tcdlcl
#ifdef ES
     $         .false., .true.,  .false., .false.,  .false., .true.,
#else
     $         .true.,  .true.,  .true.,  .true.,   .true.,  .true.,
#endif
c              preslct  preslcb  tmplct   prate     cprat    gflux
#ifdef ES
     $         .false., .false., .false., .true.,   .true.,  .true.,
#else
     $         .true.,  .true.,  .true.,  .true.,   .true.,  .true.,
#endif
c              land     icec     ugrdhag  vgrdhag   tmphag   spfhhag
     $         .true.,  .true.,  .true.,  .true.,   .true.,  .true.,
c              tmax     tmin     runof    pevpr     ugrd     vgwd
     $         .true.,  .true.,  .true.,  .true.,   .true.,  .true.,
c              hpbl     pwat     srweq    snoev    snohf     quflx
     $         .true.,  .true.,  .true.,  .true.,  .true.,   .true.,
c              qvflx    cwat     dswrftoa tcdcclm  albdo     sfcr
     $         .true.,  .true.,  .true.,  .true.,  .true.,   .true.,
c              vgtyp     veg      sltyp   alhtf    ecpy      bgrun
     $         .true.,   .true.,  .true., .true.,  .true.,   .true.,
c              psfc    aprcp   fwindx
#ifdef NOAHYDRO
     $         .true., .true., .true.
#else
     $         .true.,         .true.
#endif
#ifdef ISOTOPE
c              prate1  prate2  cprat1  cprat1  lhtfl1  lhtfl2
     $        ,.true., .true., .true., .true., .true., .true.,
c              rresv1  rresv2  pwat1   pwat2   q2m1    q2m2
     $         .true., .true., .true., .true., .true., .true.,
c              ccond   ccond1  ccond2  cond    cond1   cond2
     $         .true., .true., .true., .true., .true., .true.,
c              eqfo    eqfd    runof1  runof2  bgrun1  bgrun2
     $         .true., .true., .true., .true., .true., .true.
#endif
#endif
     $           /
c
      parameter (ipuu=33,ipuv=34,ipuo=39,ipuz=7,iput=11,ipur=52,ipua=41)
      parameter (ipuq=51,ipucl=153,ipuo3=154,ipupr=152)
c
      logical   lbm(IDIMF*JDIMF)
      dimension ipusun(nsun)
      dimension itlsun(nsun),il1sun(nsun),il2sun(nsun)
      dimension kslp(2)
      dimension iens(5)
c
      save iens
      data iens/1,0,0,1,255/
c
      logical   lppr,lpcl,lpo3

      real, allocatable :: fxs(:,:,:),fxy(:,:,:)
      real      fxsq(IDIM2,JDIM2,levh_)
c
#ifdef SMTHBZ
      real, allocatable :: rltbf(:,:)
#endif

#ifdef MP
      real, allocatable :: fullat(:,:),fullon(:,:)
#endif
#ifdef RSM
      dimension fm2(*)
#endif
c
#ifdef MP
      dimension nstr(0:npes_-1),nend(0:npes_-1)
#endif
      dimension lgrib(ntlen)
#ifdef RSM
      character grib(30+lenpds+lengds+IDIMF*JDIMF*(32+1)/8,ntlen)
#else
      dimension gout(IDIMO,JDIMO)
      character grib(30+lenpds+lengds+IDIMO*JDIMO*(32+1)/8,ntlen)
#endif
c
      dimension iens2(5)
c
      parameter(pt=30.,icen=7,icen2=0,igen=99)
#ifdef MP
c
      real, allocatable :: work(:,:,:)
#endif
c
      data ifp/0/
      save ifp
c

      character*6 fni
      character*80 fno
#ifdef MP
      character*4 cpe
#endif
      dimension ak5(levs_+1),bk5(levs_+1)
      dimension sihyb(IDIM2,levs_+1),slhyb(IDIM2,levs_)



#ifdef NISLQ
#include <comnislq.h>
#endif


c
#ifdef MP
#ifdef PGBGATHER
      integer icpp
#endif
#endif
c
      ncloud = ncldg_  
      icloud = ntrac_ + 1 

      lpo3=.false.
      lpcl=.false.
      lppr=.false.
#ifndef O3C
      lpo3=ntrac_.ge.2
#endif
#ifndef CO2PROG
      lpcl=ncldg_.ge.1
#else
      lpcl=.true.
#endif
      lppr=ncldg_.ge.2
#ifdef ISOTOPE
      lpo3=.true.
      lpcl=.true.
      lppr=.true.
#endif
c
c  set constants
c
      call idsdef(2,ids)
c
#ifdef MP
#ifndef MPIGRIB
      if(ntlen.ge.npes_) then
        nblock=max(ntlen/npes_,1)
        modnb=mod(ntlen,npes_)
        nlenmx=0
        nn=1
        do n=0,npes_-1
          nstr(n)=nn
          if(n.le.modnb) then
            nblen=nblock+1
          else
            nblen=nblock
          endif
          nend(n)=min(nn+nblen-1,ntlen)
          nlen=nend(n)-nstr(n)+1
          nlenmx=max(nlenmx,nlen)
          nn=nn+nlen
        enddo
      else
        nn=1
        do n=0,ntlen-1
          nstr(n)=nn
          nend(n)=nn
          nn=nn+1
        enddo
        do n=ntlen,npes_-1
          nstr(n)=0
          nend(n)=-1
        enddo
        nlenmx=1
      endif
      allocate (work(IDIMF,JDIMF,nlenmx))
#else
      allocate (work(IDIMF,JDIMF,1))
#endif
#endif
c
#ifdef RSM
      if(.not.lbaseout) then
        fni='r_pgb'
        nchi=5
      else
        fni='r_bpgb'
        nchi=6
      endif
#else
      fni='pgb'
      nchi=3
#endif
      call fnam(fni,nchi,fhour,fno,ncho)
#ifdef MP
#ifndef MPIGRIB
#ifndef PGBGATHER
      if(mype.lt.10) then
        ndig=2
        write(cpe,'(1h0,i1)') mype
      elseif(mype.lt.100) then
        ndig=2
        write(cpe,'(i2)') mype
      elseif(mype.lt.1000) then
        ndig=3
        write(cpe,'(i3)') mype
      elseif(mype.lt.10000) then
        ndig=4
        write(cpe,'(i4)') mype
      endif
      fno=fno(1:ncho)//'_pe'//cpe
      ncho=ncho+3+ndig
#endif
#endif
#endif
#ifdef MP
#ifdef MPIGRIB
      if(mype.eq.master) then
#endif
#ifdef PGBGATHER
      if(mype.eq.master) then
#endif
#ifndef MPIGRIB
#ifndef PGBGATHER
      if(mype.le.ntlen) then
#endif
#endif
#endif
      open(unit=n1,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      print *,' error in opening file ',fno(1:ncho)
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  901 continue 
      print *,' opening file ',fno(1:ncho)
#ifdef MP
#ifdef MPIGRIB
      endif
#endif
#ifdef PGBGATHER
      endif
#endif
#ifndef MPIGRIB
#ifndef PGBGATHER
      endif
#endif
#endif
#endif
c
c  set both input and output indices
c
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,' set both input and output indices '
#ifdef MP
      endif
#endif
      ksz=1
      ksd=1+levs_
      kst=1+2*levs_
      ksq=1+3*levs_
      kscl=1+4*levs_
      kspr=1+5*levs_
      kso3=1+6*levs_
      kspsx=1+7*levs_
      kspsy=2+7*levs_
      ksu=3+7*levs_
      ksv=3+8*levs_
      ksps=3+9*levs_
      ksgz=4+9*levs_
      ksgzx=5+9*levs_
      ksgzy=6+9*levs_
      nflds=6+9*levs_
c
      kpz=1
      kpu=1+ko_
      kpv=1+2*ko_
      kpt=1+3*ko_
      kpo=1+4*ko_
      kpr=1+5*ko_
      kpq=1+6*ko_
      kpa=1+7*ko_
      kpo3=1+8*ko_
      kpcl=1+9*ko_
      kppr=1+10*ko_
      kptu=11*ko_+1
      kptv=11*ko_+kt_+1
      kptt=11*ko_+2*kt_+1
      kptr=11*ko_+3*kt_+1
      kpsun=11*ko_+4*kt_+1
c
      allocate (fxs(IDIM2,JDIM2,nflds))
      allocate (fxy(IDIM2,JDIM2,ntlen))
c
c  set some parameters
c
      jfhour=nint(fhour)
c
      do nn=1,ntlen
        mxbit(nn)=16
        ifhr(nn)=jfhour
        ithr(nn)=0
        iavg(nn)=10
      enddo
c
      nfldp=nupa*ko_+nupt*kt_+nsun
      do k=1,ko_
        pokpa(k)=po(k)/10.
        if(float(nint(po(k))).eq.po(k).or.po(k).gt.655.) then
          ipo(k)=100
          npo(k)=nint(po(k))
        else
          ipo(k)=120
          npo(k)=nint(po(k)*100.)
        endif
      enddo
      ptkpa=pt/10.
      do k=1,kt_
        npt1(k)=k*pt
        npt2(k)=(k-1)*pt
      enddo
      kpmu=isrchflt(ko_,po,1,pot(ipuu))
      kpmv=isrchflt(ko_,po,1,pot(ipuv))
      kpmo=isrchflt(ko_,po,1,pot(ipuo))
      kpmz=isrchflt(ko_,po,1,pot(ipuz))
      kpmt=isrchflt(ko_,po,1,pot(iput))
      kpmr=isrchflt(ko_,po,1,pot(ipur))
      kpmq=isrchflt(ko_,po,1,pot(ipuq))
      kpma=isrchflt(ko_,po,1,pot(ipua))
      kpmcl=isrchflt(ko_,po,1,pot(ipucl))
      kpmpr=isrchflt(ko_,po,1,pot(ipupr))
      kpmo3=isrchflt(ko_,po,1,pot(ipuo3))
      call sunprm(ko_,po,lpcl,lppr,ipusun,itlsun,il1sun,il2sun,
     1            kslp,kli)
c
c  set some grib parameters
c
      do nn=1,nfldp
        ipu(nn)=0
      enddo
      do nn=kpu,kpu+kpmu-1
        ipu(nn)=ipuu
      enddo
      do nn=kpv,kpv+kpmv-1
        ipu(nn)=ipuv
      enddo
      do nn=kpo,kpo+kpmo-1
        ipu(nn)=ipuo
      enddo
      do nn=kpz,kpz+kpmz-1
        ipu(nn)=ipuz
      enddo
      do nn=kpt,kpt+kpmt-1
        ipu(nn)=iput
      enddo
      do nn=kpr,kpr+kpmr-1
        ipu(nn)=ipur
      enddo
      do nn=kpq,kpq+kpmq-1
        ipu(nn)=ipuq
      enddo
      do nn=kpa,kpa+kpma-1
        ipu(nn)=ipua
      enddo
      do nn=kpcl,kpcl+kpmcl-1
        if(lpcl) then
          ipu(nn)=ipucl
        else
          ipu(nn)=0
        endif
      enddo
      do nn=kppr,kppr+kpmpr-1
        if(lppr) then
          ipu(nn)=ipupr
        else
          ipu(nn)=0
        endif
      enddo
      do nn=kpo3,kpo3+kpmo3-1
        if(lpo3) then
          ipu(nn)=ipuo3
        else
          ipu(nn)=0
        endif
      enddo
      do nn=kptu,kptu+kt_-1
        ipu(nn)=ipuu
      enddo
      do nn=kptv,kptv+kt_-1
        ipu(nn)=ipuv
      enddo
      do nn=kptt,kptt+kt_-1
        ipu(nn)=iput
      enddo
      do nn=kptr,kptr+kt_-1
        ipu(nn)=ipur
      enddo
      do nn=kpsun,kpsun+nsun-1
        ipu(nn)=ipusun(nn-kpsun+1)
      enddo
c
      do nn=kpu,kpu+ko_-1
        itl(nn)=ipo(nn-kpu+1)
      enddo
      do nn=kpv,kpv+ko_-1
        itl(nn)=ipo(nn-kpv+1)
      enddo
      do nn=kpo,kpo+ko_-1
        itl(nn)=ipo(nn-kpo+1)
      enddo
      do nn=kpz,kpz+ko_-1
        itl(nn)=ipo(nn-kpz+1)
      enddo
      do nn=kpt,kpt+ko_-1
        itl(nn)=ipo(nn-kpt+1)
      enddo
      do nn=kpr,kpr+ko_-1
        itl(nn)=ipo(nn-kpr+1)
      enddo
      do nn=kpq,kpq+ko_-1
        itl(nn)=ipo(nn-kpq+1)
      enddo
      do nn=kpa,kpa+ko_-1
        itl(nn)=ipo(nn-kpa+1)
      enddo
      do nn=kpo3,kpo3+ko_-1
        itl(nn)=ipo(nn-kpo3+1)
      enddo
      do nn=kpcl,kpcl+ko_-1
        itl(nn)=ipo(nn-kpcl+1)
      enddo
      do nn=kppr,kppr+ko_-1
        itl(nn)=ipo(nn-kppr+1)
      enddo
      do nn=nupa*ko_+1,nupa*ko_+nupt*kt_
        itl(nn)=116
      enddo
      do nn=kpsun,kpsun+nsun-1
        itl(nn)=itlsun(nn-kpsun+1)
      enddo
      do nn=1,nupa*ko_+nupt*kt_+nsun
        il1(nn)=0
      enddo
      do nn=kptu,kptu+kt_-1
        il1(nn)=npt1(nn-kptu+1)
      enddo
      do nn=kptv,kptv+kt_-1
        il1(nn)=npt1(nn-kptv+1)
      enddo
      do nn=kptt,kptt+kt_-1
        il1(nn)=npt1(nn-kptt+1)
      enddo
      do nn=kptr,kptr+kt_-1
        il1(nn)=npt1(nn-kptr+1)
      enddo
      do nn=kpsun,kpsun+nsun-1
        il1(nn)=il1sun(nn-kpsun+1)
      enddo
      do nn=kpu,kpu+ko_-1
        il2(nn)=npo(nn-kpu+1)
      enddo
      do nn=kpv,kpv+ko_-1
        il2(nn)=npo(nn-kpv+1)
      enddo
      do nn=kpo,kpo+ko_-1
        il2(nn)=npo(nn-kpo+1)
      enddo
      do nn=kpz,kpz+ko_-1
        il2(nn)=npo(nn-kpz+1)
      enddo
      do nn=kpt,kpt+ko_-1
        il2(nn)=npo(nn-kpt+1)
      enddo
      do nn=kpr,kpr+ko_-1
        il2(nn)=npo(nn-kpr+1)
      enddo
      do nn=kpq,kpq+ko_-1
        il2(nn)=npo(nn-kpq+1)
      enddo
      do nn=kpa,kpa+ko_-1
        il2(nn)=npo(nn-kpa+1)
      enddo
      do nn=kpo3,kpo3+ko_-1
        il2(nn)=npo(nn-kpo3+1)
      enddo
      do nn=kpcl,kpcl+ko_-1
        il2(nn)=npo(nn-kpcl+1)
      enddo
      do nn=kppr,kppr+ko_-1
        il2(nn)=npo(nn-kppr+1)
      enddo
      do nn=kptu,kptu+kt_-1
        il2(nn)=npt2(nn-kptu+1)
      enddo
      do nn=kptv,kptv+kt_-1
        il2(nn)=npt2(nn-kptv+1)
      enddo
      do nn=kptt,kptt+kt_-1
        il2(nn)=npt2(nn-kptt+1)
      enddo
      do nn=kptr,kptr+kt_-1
        il2(nn)=npt2(nn-kptr+1)
      enddo
      do nn=kpsun,kpsun+nsun-1
        il2(nn)=il2sun(nn-kpsun+1)
      enddo
c
      if(.not.lskipgb) then
c
c  convert wave to grid
c
#ifdef RSM
#define COF2GRD rcof2grd
#else
#define COF2GRD cof2grd
#endif
#ifdef RSM
      call COF2GRD(fgz,fq,fte,fuu,fvv,frq,fm2,
#else
      call COF2GRD(fgz,fq,fte,fuu,fvv,frq,
#endif
     $             fxs(1,1,ksgz),fxs(1,1,ksgzx),fxs(1,1,ksgzy),
     $             fxs(1,1,ksps),fxs(1,1,kspsx),fxs(1,1,kspsy),
     $             fxs(1,1,kst ),
     $             fxs(1,1,ksu ),fxs(1,1,ksv  ),
     $             fxs(1,1,ksd ),fxs(1,1,ksz  ),
     $             fxsq(1,1,1 )
#ifdef RSM
     $            ,lbaseout)
#else
     $             )
#endif
c

#ifdef NISLQ
        do k = 1, levh_
          do j = 1, JDIM2
            do i = 1, IDIM2
              fxsq(i,j,k) = slq_q2(i,k,j)
            enddo
          enddo
        enddo
#endif


      do k=1,levs_
        do j=1,JDIM2
          do i=1,IDIM2
            fxs(i,j,ksq+k-1)=fxsq(i,j,k)
          end do
        end do
      end do

      if(lpcl) then
#ifndef CO2PROG
#ifdef ISOTOPE
        kc=levs_*2+1      ! 18o (tracer 3)
#else
        kc=ntrac_*levs_+1      ! cloud water
#endif
#else
        kc=levs_*2+1      ! co2 (tracer 3)
#endif
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
#ifdef ISOTOPE               
              fxs(i,j,kscl+k-1)=fxsq(i,j,kc+k-1)
#else
              fxs(i,j,kscl+k-1)=max(0.0,fxsq(i,j,kc+k-1))
#endif
            end do
          end do
        end do        
#if ( _ncldg_ >= 4 )
        kc=(ntrac_+1)*levs_+1  ! cloud ice
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
              fxs(i,j,kscl+k-1)= fxs(i,j,kscl+k-1)
     &           + max(0.0,fxsq(i,j,kc+k-1))
            end do
          end do
        end do        
#endif
      end if
      if(lppr) then
#if ( _ncldg_ >= 4 )
        kc=(ntrac_+2)*levs_+1   ! rain water
#else
        kc=(ntrac_+1)*levs_+1   ! rain water
#endif
#ifdef ISOTOPE
        kc=levs_*3+1   ! hdo (tracer 4)
#endif
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
#ifdef ISOTOPE
              fxs(i,j,kspr+k-1)=fxsq(i,j,kc+k-1)
#else
              fxs(i,j,kspr+k-1)=max(0.0,fxsq(i,j,kc+k-1))
#endif
            end do
          end do
        end do        
#if ( _ncldg_ >= 4 )
        kc=(ntrac_+3)*levs_+1  ! snow
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
              fxs(i,j,kspr+k-1)= fxs(i,j,kspr+k-1)
     &           + max(0.0,fxsq(i,j,kc+k-1))
            end do
          end do
        end do        
#endif
#if ( _ncldg_ >= 5 )
        kc=(ntrac_+4)*levs_+1  ! graupel
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
              fxs(i,j,kspr+k-1)= fxs(i,j,kspr+k-1)
     &           + max(0.0,fxsq(i,j,kc+k-1))
            end do
          end do
        end do        
#endif
      end if
      if(lpo3) then
        kc=levs_+1  ! o3 is supposed as second fxsq
        do k=1,levs_
          do j=1,JDIM2
            do i=1,IDIM2
              fxs(i,j,kso3+k-1)=fxsq(i,j,kc+k-1)
            end do
          end do
        end do        
      end if

c
c  convert virtual temp to real temp
c
      do k=1,levs_
#ifndef MP
        do j=1,JDIM2
          do i=1,IDIM2
            w2 = 1.+ fvirt*fxsq(i,j,k)
            if(ncloud.gt.0) then
              w1 = 0.
              do ic = icloud,ntotal_
                kc = (ic-1)*levs_ + k
                w1 = max(fxsq(i,j,kc),qmin) + w1
              enddo
              w2 = w2 - w1
            endif

            fxs(i,j,kst+k-1)=fxs(i,j,kst+k-1)/w2

          enddo
        enddo
#else
        do j=1,latlen(mype)
          do i=1,lonlen(mype)*2
            fxs(i,j,kst+k-1)=fxs(i,j,kst+k-1)/
     1                      (1.+fvirt*fxs(i,j,ksq+k-1))
          enddo
        enddo
#endif
      enddo
c
c  convert ln of ps to real ps
c
#ifndef MP
      do j=1,JDIM2
        do i=1,IDIM2
          fxs(i,j,ksps)=exp(fxs(i,j,ksps))
        enddo
      enddo
#else
      do j=1,latlen(mype)
        do i=1,lonlen(mype)*2
          fxs(i,j,ksps)=exp(fxs(i,j,ksps))
        enddo
      enddo
#endif
#ifdef DBG
      call maxmin(fxs(1,1,ksgz ),IDIM2*JDIM2, 1,1, 1,'gz')
      call maxmin(fxs(1,1,ksgzx),IDIM2*JDIM2, 1,1, 1,'dgz/dx')
      call maxmin(fxs(1,1,ksgzy),IDIM2*JDIM2, 1,1, 1,'dgz/dy')
      call maxmin(fxs(1,1,ksps ),IDIM2*JDIM2, 1,1, 1,'ps')
      call maxmin(fxs(1,1,kspsx),IDIM2*JDIM2, 1,1, 1,'dps/dx')
      call maxmin(fxs(1,1,kspsy),IDIM2*JDIM2, 1,1, 1,'dps/dy')
      call maxmin(fxs(1,1,kst  ),IDIM2*JDIM2,levs_,1,levs_,'tmp')
      call maxmin(fxs(1,1,ksu  ),IDIM2*JDIM2,levs_,1,levs_,'u')
      call maxmin(fxs(1,1,ksv  ),IDIM2*JDIM2,levs_,1,levs_,'v')
      call maxmin(fxs(1,1,ksd  ),IDIM2*JDIM2,levs_,1,levs_,'div')
      call maxmin(fxs(1,1,ksz  ),IDIM2*JDIM2,levs_,1,levs_,'vot')
      call maxmin(fxs(1,1,ksq  ),IDIM2*JDIM2,levs_,1,levs_,'q')
      if(lpcl)
     * call maxmin(fxs(1,1,kscl ),IDIM2*JDIM2,levs_,1,levs_,'cl')
      if(lppr)
     * call maxmin(fxs(1,1,kspr ),IDIM2*JDIM2,levs_,1,levs_,'pr')
      if(lpo3)
     * call maxmin(fxs(1,1,kso3 ),IDIM2*JDIM2,levs_,1,levs_,'o3')
#endif
c
c  loop over groups of latitudes
c
c  compute auxiliary quantities on sigma and interpolate to pressure
c  and compute sundry fields and pack for transpose in parallel
c
#ifdef MP
c
c  fill undefined partial array points
c
        do n=1,nflds
          do j=1,JDIM2
            if(j.gt.latlen(mype)) then
              do i=1,IDIM2
                fxs(i,j,n)=fxs(1,1,n)
              enddo
            endif
          enddo
          do i=1,IDIM2
            if(i.gt.lonlen(mype)*2) then
              do j=1,JDIM2
                fxs(i,j,n)=fxs(1,1,n)
              enddo
            endif
          enddo
        enddo
#ifdef RSM
        do j=1,JDIM2
          if(j.gt.latlen(mype)) then
            do i=1,IDIM2
              flat(i,j)=flat(1,1)
            enddo
          endif
        enddo
        do i=1,IDIM2
          if(i.gt.lonlen(mype)*2) then
            do j=1,JDIM2
              flat(i,j)=flat(1,1)
            enddo
          endif
        enddo
#endif
#endif
c
        j1=1
        j2=JDIM2
        ijlen=IDIM2*JDIM2
c
        do j=j1,j2
          jj=j-j1+1
#ifdef HYBRID          
          do i=1,IDIM2
            do k=1,levs_+1
              sihyb(i,k)=ak5(levs_-k+2)/fxs(i,jj,ksps)+bk5(levs_-k+2)
            enddo
            do k=1,levs_
              slhyb(i,k)=0.5*(sihyb(i,k)+sihyb(i,k+1))
            enddo
          enddo
#else
          do i=1,IDIM2
            do k=1,levs_+1
              sihyb(i,k)=si(k)
            enddo
            do k=1,levs_
              slhyb(i,k)=sl(k)
            enddo
          enddo          
#endif
          call getrh(IDIM2,ijlen,levs_,slhyb,
     &             fxs(1,jj,ksps),fxs(1,jj,ksq),fxs(1,jj,kst),
     &             qsxs,rxs)
          call omega(IDIM2,ijlen,levs_,sihyb,slhyb,
     &             fxs(1,jj,ksps),fxs(1,jj,kspsx),fxs(1,jj,kspsy),
     &             fxs(1,jj,ksd ),fxs(1,jj,ksu  ),fxs(1,jj,ksv  ),
     &             oxs,osxs)
          call hydro(IDIM2,ijlen,levs_,sihyb,slhyb,
     &             fxs(1,jj,ksgz),fxs(1,jj,kst),fxs(1,jj,ksq),
     &             zxs,zxi)
          call sig2p(IDIM2,ijlen,levs_,sihyb,slhyb,fxs(1,jj,ksps),
     &             fxs(1,jj,ksu),fxs(1,jj,ksv),oxs,
     &             zxs,zxi,fxs(1,jj,kst),rxs,fxs(1,jj,ksq),
     &             fxs(1,jj,ksz),
     &             fxs(1,jj,kscl),fxs(1,jj,kspr),fxs(1,jj,kso3),
     &             ko_,pokpa,lpcl,lppr,lpo3,
     &             fxp(1,kpu),fxp(1,kpv),fxp(1,kpo),
     &             fxp(1,kpz),fxp(1,kpt),fxp(1,kpr),
     &             fxp(1,kpq),fxp(1,kpa),fxp(1,kpcl),fxp(1,kppr),
     &             fxp(1,kpo3))
          call sig2pt(IDIM2,ijlen,levs_,sihyb,slhyb,fxs(1,jj,ksps),
     &              fxs(1,jj,ksu),fxs(1,jj,ksv),
     &              fxs(1,jj,kst),fxs(1,jj,ksq),qsxs,
     &              kt_,ptkpa,
     &              fxp(1,kptu),fxp(1,kptv),
     &              fxp(1,kptt),qxp,fxp(1,kptr))
          call sundry(IDIM2,ijlen,levs_,kslp,kli,lpcl,lppr,
     &              flat(1,j),sihyb,slhyb,kt_,ptkpa,
     &              fxs(1,jj,ksgz),fxs(1,jj,ksps),osxs,
     &              fxs(1,jj,ksu),fxs(1,jj,ksv),oxs,
     &              fxs(1,jj,kst),rxs,fxs(1,jj,ksq),qsxs,
     &              fxp(1,kptt),qxp,
     &              fxp(1,kpz),fxp(1,kpt),
     &              fxs(1,jj,kscl),fxs(1,jj,kspr),
     &              fxp(1,kpsun))
c
c compute absolute vorticity
c
          aomega=2.0*acos(-1.0)/(24.*60.*60.)
          do nn=1,nfldp
            do i=1,IDIM2 
              fxy(i,jj,nn)=fxp(i,nn) 
            enddo
          enddo
          do k=1,ko_
            do i=1,IDIM2
               fxy(i,jj,kpa-1+k)= fxy(i,jj,kpa-1+k)
#ifdef RSM
     1                     +2.0*aomega*sin(flat(i,j))
#else
     1                     +2.0*aomega*sin(flat(1,j))
#endif
             enddo
          enddo
c
#ifdef COUPLE_RSM_AND_ROMS
          do i=1,IDIM2
             romssgz(i,jj)= fxs(i,jj,ksgz)
          enddo
#endif
c
        enddo
c
c  endif for if(.not.lskipgb)
c
      endif 
c
#ifdef RSM
      n=ntlen-nflx
      call addflx(fxy,n,mxbit,ipu,itl,il1,il2,ifhr,ithr,iavg)
      if(n.ne.ntlen) then
          print *,'wrong number of flux records'
          call abort
      endif
#ifdef SMTHBZ
      allocate(rltbf(igrd1_,jgrd1_))
#ifdef MP
      call rmpgp2f(rltb,igrd12p_,jgrd12p_,rltbf,igrd12_,jgrd12_,1)
      call rmpbcastr(rltbf,igrd12_*jgrd12_)
#else
      do ij=1,igrd12_*jgrd12_
        rltbf(ij,1)=rltb(ij)
      enddo
#endif
      call shalfo(rltbf,1)
#endif
#endif
c
#ifdef RSM
      proj=rproj
      delx=rdelx
      dely=rdely
      rlat1r=rlat1*rad
      rlon1r=rlon1*rad
      rlat2r=rlat2*rad
      rlon2r=rlon2*rad
      if( proj.eq. 0. ) then
        idrt=1                    ! mercater
        ortru=rtruth
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,' mercater projection.'
#ifdef MP
        endif
#endif
      elseif( abs(proj).eq.1. ) then
        idrt=5                    ! polar projection
        ortru=rorient
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,' polar projection.'
#ifdef MP
        endif
#endif
      else
        idrt=0
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,' undefined map projection.'
#ifdef MP
        endif
#endif
      endif
#else
      proj=0.
      delx=0.
      dely=0.
      ortru=0.
      rlat1r=0.
      rlon1r=0.
      rlat2r=0.
      rlon2r=0.
#endif
#ifndef MP
c
c  unpack transposed fields and interpolate to output grid
c  and round to the number of bits and engrib the field in parallel
c
      k1=1
      k2=ntlen
#ifdef ORIGIN
c$doacross share(k1,k2,pot,fxy,
c$&       ipu,itl,il1,il2,icen,igen,idate,lenpds,jfhour,ids,mxbit,
c$&       grib,lgrib,itruth,idelx,idely,
c$&       ilat1,ilon1,ilat2,ilon2,idrt,
c$&       rmax,rmin,icen2,rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj),
c$& local(k,kan,ierr,pok,poktop,i,j)
#endif
#ifdef CRAY
cmic$ do all
cmic$1 shared(k1,k2,pot,igrd1,fxy)
cmic$1 shared(ipu,itl,il1,il2,icen,igen,idate,lenpds,jfhour,ids,mxbit)
cmic$1 shared(grib,lgrib,itruth,idelx,idely)
cmic$1 shared(ilat1,ilon1,ilat2,ilon2,idrt)
cmic$1 shared(rmax,rmin,icen2,rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj)
cmic$1 private(k,kan,ierr,pok,poktop,i,j)
#endif
#ifdef OPENMP
!$omp parallel do private(k,kan,ierr,pok,poktop,i,j)
#endif
#endif /* ifndef MP */
#ifdef MP
#define FXY work
#else
#define FXY fxy
#endif
c
#ifndef MP
        mype=0
        do k=k1,k2
          kan=k-k1+1
#else
#ifdef MPIGRIB
        do k=1,ntlen
          kan=1
#ifndef RSM
          call mpgp2f(fxy(1,1,k),lonf2p_,latg2p_,work,lonf2_,latg2_,1)
#else
          call MPSYNALL
          call rmpgp2f(fxy(1,1,k),igrd12p_,jgrd12p_,work,
     1                 igrd12_,jgrd12_,1)
#endif
          if(mype.eq.master) then
#else  /* ifndef MPIGRIB */
#ifdef RSM
#define MPGPFK2FPK rmpgpfk2fpk
#else
#define MPGPFK2FPK mpgpfk2fpk
#endif
          call MPSYNALL
          call MPGPFK2FPK(fxy,IDIM2,JDIM2,ntlen,
     1                    work,IDIMF2,JDIMF2,nlenmx,
     2                    nstr,nend)
          do k=nstr(mype),nend(mype)
            kan=k-nstr(mype)+1
#endif  /* ifndef MPIGRIB */
c
#endif  /* ifdef MP */
            lsel=nsel(k)
            if(lskipgb.and.k.le.npgblen) lsel=.false.
            lgrib(kan)=0
            if(ipu(k).le.0.or.itl(k).eq.107.or.itl(k).eq.116) then
              lsel=.false.
            endif
            if(lsel) then
#ifdef DBG
#ifdef CRAY
cmic$ guard
#endif
              rmax=FXY(1,1,kan)
              rmin=FXY(1,1,kan)
#ifdef MP
              do j=1,JDIMF2
#ifdef REDUCE_GRID
                lonsd2=lonfdp(j,mype)*2
#else
                lonsd2=IDIMF2
#endif
                do i=1,lonsd2
#else
              do j=1,JDIM2
#ifdef REDUCE_GRID
                lonsd2=lonfd(latdef(j))*2
#else
                lonsd2=IDIM2
#endif
                do i=1,lonsd2
#endif  /* #else ifdef MP  */
                  rmax=max(FXY(i,j,kan),rmax)
                  rmin=min(FXY(i,j,kan),rmin)
                enddo
              enddo
#endif
#ifndef RSM
              call rowsep(FXY(1,1,kan))
              call gau2ll(FXY(1,1,kan),IDIMF,JDIMF,
     &                    0.,90.,360./float(IDIMO),180./float(JDIMO-1),
     &                    gout,IDIMO,JDIMO)
              call gribit(gout        ,lbm,   0,io_,jo_,mxbit(k),90.,
#else
              call shalfo(FXY(1,1,kan),1)
#ifdef SMTHBZ
              nitr=15
              do itr=1,nitr
                do j=2,jgrd1_-1
                  do i=2,igrd1_-1
                    if(rltbf(i,j).ne.1.) then
                      FXY(i,j,kan)=(FXY(i-1,j,kan)+FXY(i,j-1,kan)+
     &                            FXY(i+1,j,kan)+FXY(i,j+1,kan)+
     &                            FXY(i,j,kan)*4.0)*0.125
                    endif
                  enddo
                enddo
                do i=1,igrd1_
                  FXY(i,     1,kan)=FXY(i,       2,kan)
                  FXY(i,jgrd1_,kan)=FXY(i,jgrd1_-1,kan)
                enddo
                do j=1,jgrd1_
                  FXY(1     ,j,kan)=FXY(2       ,j,kan)
                  FXY(igrd1_,j,kan)=FXY(igrd1_-1,j,kan)
                enddo
              enddo
#endif
              call gribit(FXY(1,1,kan),lbm,idrt,igrd1_,jgrd1_,
     &                  mxbit(k),0.,
#endif
     &                  lenpds,2,icen,igen,0,
     &                  ipu(k),itl(k),il1(k),il2(k),
     &                  idate(4),idate(2),idate(3),idate(1),
     &                  1,ifhr(k),ithr(k),iavg(k),
     &                  0,0,icen2,ids(ipu(k)),iens,
     &                  rlat1r,rlon1r,rlat2r,rlon2r,delx,dely,
     '                  ortru,proj,grib(1,kan),lgrib(kan),ierr)
#ifdef DBG
              write(*,'("mype=",i3," k=",i4," max=",
     &              e17.8," min=",e17.8," lgrib=",i8)') 
     &              mype,k,rmax,rmin,lgrib(kan)
#endif
#ifdef MPIGRIB 
              if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kpcl) lsel=.false.
              if((.not.lpcl).and.k.ge.kpcl.and.k.lt.kppr) lsel=.false.
              if((.not.lppr).and.k.ge.kppr.and.k.lt.kptu) lsel=.false.
              if(lgrib(kan).gt.0.and.lsel) then
                call wryte(n1,lgrib(kan),grib(1,kan))
#ifdef DBG
                print *,' grib1 written to ',n1,'length=',lgrib(kan)
#endif
              endif
#endif  /* ifdef MPIGRIB */
#ifdef CRAY
cmic$ guard
#endif
            endif
        enddo
#ifdef MP
#ifdef MPIGRIB
            endif
#endif  /* ifdef MPIGRIB  */
#endif  /* ifdef MP */
c
#ifdef MP
#ifdef PGBGATHER
          icpp=30+lenpds+lengds+IDIMO*JDIMO*(32+1)/8
          call mpgfpk2fk(grib,icpp,
     &                   lgrib,nstr,nend,grib,ntlen)
c
#endif
#endif
#ifndef MPIGRIB
#ifdef MP 
#ifndef PGBGATHER
        do k=nstr(mype),nend(mype)
          kan=k-nstr(mype)+1
#else
        if(mype.eq.master) then
        do k=1,ntlen
          kan=k
#endif
#else
        do k=k1,k2
          kan=k-k1+1
#endif
          lsel=nsel(k)
          if(lskipgb.and.k.le.npgblen) lsel=.false.
#ifndef ISOTOPE
          if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kpcl) lsel=.false.
#endif
          if((.not.lpcl).and.k.ge.kpcl.and.k.lt.kppr) lsel=.false.
          if((.not.lppr).and.k.ge.kppr.and.k.lt.kptu) lsel=.false.
          if(lsel) then
            call wryte(n1,lgrib(kan),grib(1,kan))
#ifdef DBG
            print *,' grib1 written to ',n1,' of length ',lgrib(kan)
#endif
          endif
        enddo
#endif  /* ifndef MPIGRIB */
#ifdef MP
#ifdef PGBGATHER
        endif
#endif
#endif
cifndef MP
c     enddo
cendif
      deallocate (fxs,fxy)
c
#ifdef MP
#ifndef MPIGRIB
      deallocate (work)
#endif
#ifdef SMTHBZ
#ifdef RSM
      deallocate (rltbf)
#endif
#endif
#endif
c
c  for global model, separate pgb and flx files
c
#ifndef RSM
#ifdef MP
#ifdef MPIGRIB
      if(mype.eq.master) then
#endif
#ifdef PGBGATHER
      if(mype.eq.master) then
#endif
#endif
      close (n1)
#ifdef MP
#ifdef MPIGRIB
      endif
#endif
#ifdef PGBGATHER
      endif
#endif
#endif
#ifndef COUPLE_OUT
      call wriflx(n1)
#else
      call wriflx(n1,nnp)
#endif
#else  /* ifdef RSM */
      close(n1)

#endif
c
#endif   /* ifdef MRG_POST */
      return

      end
