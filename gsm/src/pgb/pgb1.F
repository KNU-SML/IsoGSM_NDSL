#define IO io
#define JCAP jcap
#define JO jo
#define KO ko
#define LEVS levs
#define NCLDG ncldg
#define NCPUS ncpus
#define NTOTAL ntotal
#define NTRAC ntrac
#define RD rd
#define RV rv
#define CP cp
#include <define.h>
      subroutine pgb1(fhour,idate,nsig,si,sl,po,
     &            npgb,ncpus,ids,pot,icen,icen2,igen,ntrace,ncloud)
#include <paramodel.h>
c-paramodel
c
c  in this routine wave number is an input model resoultion
c
      parameter(io=io_, jo=jo_, ko=ko_, mxbit=16)
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)
      parameter(jcap=jcap_,levs=levs_)
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)
      parameter(nflds=5*levs_+ntotal_*levs_+6)
c
c$$$  subprogram documentation block
c
c subprogram: pgb1           transforms sigma to pressure grib
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: transforms a sigma spectral file to pressure grib1.
c   one latitude slice at a time, first sigma grid data
c   is transformed from sigma spectral coefficients.
c   the input data consists of vorticity, divergence, wind components,
c   temperature and specific humidity on the sigma surfaces as well as
c   surface pressure and orography and their horizontal gradients,
c   plus the potential of the bottom pressure gradient force.
c   relative humidity, vertical velocity and geopotential heights
c   are computed on the sigma surfaces and then interpolated to pressure
c   along with wind and temperature.  sundry fields are also computed.
c   the output data is quarterpacked and transposed to horizontal fields
c   which are then interpolated to the output grid and rounded
c   and packed into grib messages and written to the pressure grib file.
c
c program history log:
c   92-10-31  iredell
c
c usage:    call pgb1(fhour,idate,nsig,si,sl,
c    &                npgb,ncpus,ids,pot,icen,icen2,igen)
c   input arguments:
c     fhour        real forecast hour
c     idate        integer (4) date
c     nsig         integer unit from which to read sigma file
c     si           real (levs+1) sigma interface values
c     sl           real (levs) sigma full level values
c     po           real (ko) mandatory pressures in kpa
c     npgb         integer unit to which to write grib messages
c     ncpus        integer number of cpus over which to distribute work
c     ids          integer (255) decimal scaling
c     pot          real (255) topmost pressure in kpa
c     icen         integer forecast center identifier
c     icen2        integer forecast subcenter identifier
c     igen         integer generating model identifier
c
c subprograms called:
c   rdss         read sigma coefficients
c   sunprm       set parameters for sundry fields
c   mpfdef       set defaults for field parameter identifiers
c   trss         transform sigma coefficients
c   getrh        compute relative humidity
c   omega        compute vertical velocity
c   hydro        compute geopotential heights
c   sig2p        interpolate sigma to pressure
c   sundry       compute sundry fields
c   ptranw       quarterpack and transpose data
c   ptranr       unpack quarterpacked transposed data
c   polext       extrapolate pole values
c   rowsep       separate hemispheric rows on grid
c   gribit       create grib message
c   wryte        write data by bytes
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension idate(4),si(levs+1),sl(levs),po(ko)
      dimension sihyb(io22,levs+1),slhyb(io22,levs)
      dimension ids(255),pot(255)
c
      parameter(nupa=11,nsun=4)
      parameter(lenpds=28,lengds=32)
c
      dimension ifax(20)
      dimension trig(io2),eps(nc/2),epstop(nctop/2)
      dimension ss(nc,nflds),sstop(nctop,nflds)
      dimension ipo(ko),npo(ko),pokpa(ko)
      dimension clat(johf),slat(johf),wlat(johf)
      dimension fxs(io22,nflds)
      dimension oxs(io22,levs),rxs(io22,levs)
      dimension zxs(io22,levs),zxi(io22,levs)
      dimension fxp(io22,nupa*ko+nsun)
      dimension fxy(io2,johf,nupa*ko+nsun)
c
#define DEFAULT
#ifdef DYNAMIC_ALLOC
#undef DEFAULT
      dimension lgrib(ncpus)
      character grib(30+lenpds+lengds+io*jo*(mxbit+1)/8,ncpus)
#endif
#ifdef DEFAULT
      dimension lgrib(1)
      character grib(30+lenpds+lengds+io*jo*(mxbit+1)/8,1)
#endif
      dimension ipu(nupa*ko+nsun),itl(nupa*ko+nsun)
      dimension il1(nupa*ko+nsun),il2(nupa*ko+nsun)
      dimension mpf(255)
      parameter(ipuu=33,ipuv=34,ipuo=39,ipuz=7,iput=11,ipur=52,ipus=175)
      parameter(ipuq=51,ipucl=153,ipuo3=154,ipupr=152)
      parameter(ipua=41)
      dimension ipusun(nsun)
      dimension itlsun(nsun),il1sun(nsun),il2sun(nsun)
      dimension ip1sun(nsun),ip2sun(nsun)
      dimension kslp(2)
c
      dimension iens(5)
c
      dimension ipudef(nsun),itldef(nsun)
      data ipudef/001,054,007,001/
      data itldef/001,200,001,102/
      logical  lppr,lpcl,lpo3
c
      character*3 fni
      parameter(fni='pgb',nchi=3)
      character*80 fno
#ifdef ASSIGN
      character*80 asgnstr
#endif
c
      call fnam(fni,nchi,fhour,fno,ncho)
#ifdef ASSIGN
      write(asgnstr,'(23hassign -s unblocked  u:,I2,)') npgb
      call assign('assign -R')
      call assign(asgnstr)
#endif
      open(unit=npgb,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file ',fno(1:ncho)
      call abort
  901 continue
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',npgb
c
c  set some parameters
c
      call rdss(nsig,nflds,ntrace,ncloud,nc,nctop,jcap,
     1          sl,clat,slat,wlat,trig,ifax,eps,epstop,ss,sstop)
c
      lpcl=ncloud.ge.1
      lpo3=ntrace.ge.2
      lppr=ncloud.ge.2
      jfhour=nint(fhour)
      nfldp=nupa*ko+nsun
      do k=1,ko
        pokpa(k)=po(k)/10.
#ifdef SCM
          ipo(k)=100
          npo(k)=nint(po(k))
#else
        if(float(nint(po(k))).eq.po(k).or.po(k).gt.655.) then
          ipo(k)=100
          npo(k)=nint(po(k))
        else
          ipo(k)=120
          npo(k)=nint(po(k)*100.)
        endif
#endif
      enddo
      kpmo=ko
      kpmr=ko
      kpmu=ko
      kpmv=ko
      kpmt=ko
      kpmz=ko
      kpma=ko
      kpmq=ko
      kpmcl=ko
      kpmpr=ko
      kpmo3=ko
      call sunprm(levs,ko,po,nsun,
     &            ipusun,itlsun,ip1sun,ip2sun,kslp)
      do n=1,nsun
        ipusun(n)=ipudef(n)
        itlsun(n)=itldef(n)
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  set both input and output indices
      kszs=1
      kszsx=2
      kszsy=3
      ksps=4
      kspsx=5
      kspsy=6
      ksd=7
      ksz=levs+7
      ksu=2*levs+7
      ksv=3*levs+7
      kst=4*levs+7
      ksq=5*levs+7
      kso3=ksq
      kscl=ksq
      kspr=ksq
      if(lpo3) then
         kso3=6*levs+7
         kscl=7*levs+7
         kspr=8*levs+7
      elseif(lpcl) then
         kscl=6*levs+7
         kspr=7*levs+7
      else
         kspr=6*levs+7
      endif
      kpz=1
      kpu=ko+1
      kpv=2*ko+1
      kpr=3*ko+1
      kpq=4*ko+1
      kpt=5*ko+1
      kpo=6*ko+1
      kpa=7*ko+1
      kpo3=8*ko+1
      kpcl=9*ko+1
      kppr=10*ko+1
      kpsun=11*ko+1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  set some grib parameters
      do i=1,nupa*ko
      ipu(i)=1
      enddo
      do i=kpu,kpu+kpmu-1
      ipu(i)=ipuu
      enddo
      do i=kpv,kpv+kpmv-1
      ipu(i)=ipuv
      enddo
      do i=kpo,kpo+kpmo-1
      ipu(i)=ipuo
      enddo
      do i=kpz,kpz+kpmz-1
      ipu(i)=ipuz
      enddo
      do i=kpt,kpt+kpmt-1
      ipu(i)=iput
      enddo
      do i=kpa,kpa+kpma-1
      ipu(i)=ipua
      enddo
      do i=kpr,kpr+kpmr-1
      ipu(i)=ipur
      enddo
      do i=kpq,kpq+kpmq-1
      ipu(i)=ipuq
      enddo
      do i=kpcl,kpcl+kpmcl-1
      ipu(i)=ipucl
      enddo
      do i=kppr,kppr+kpmpr-1
      ipu(i)=ipupr
      enddo
      do i=kpo3,kpo3+kpmo3-1
      ipu(i)=ipuo3
      enddo
      do i=kpsun,kpsun+nsun-1
      ipu(i)=ipusun(i-kpsun+1)
      enddo
c
      do i=kpu,kpu+ko-1
      itl(i)=ipo(i-kpu+1)
      enddo
      do i=kpv,kpv+ko-1
      itl(i)=ipo(i-kpv+1)
      enddo
      do i=kpo,kpo+ko-1
      itl(i)=ipo(i-kpo+1)
      enddo
      do i=kpz,kpz+ko-1
      itl(i)=ipo(i-kpz+1)
      enddo
      do i=kpt,kpt+ko-1
      itl(i)=ipo(i-kpt+1)
      enddo
      do i=kpa,kpa+ko-1
      itl(i)=ipo(i-kpa+1)
      enddo
      do i=kpr,kpr+ko-1
      itl(i)=ipo(i-kpr+1)
      enddo
      do i=kpq,kpq+ko-1
      itl(i)=ipo(i-kpq+1)
      enddo
      do i=kpo3,kpo3+ko-1
      itl(i)=ipo(i-kpo3+1)
      enddo
      do i=kpcl,kpcl+ko-1
      itl(i)=ipo(i-kpcl+1)
      enddo
      do i=kppr,kppr+ko-1
      itl(i)=ipo(i-kppr+1)
      enddo
      do i=kpsun,kpsun+nsun-1
      itl(i)=itlsun(i-kpsun+1)
      enddo
c
      do i=1,nupa*ko+nsun
      il1(i)=0
      enddo
      do i=kpsun,kpsun+nsun-1
c     il1(i)=il1sun(i-kpsun+1)
      il1(i)=0
      enddo
c
      do i=kpu,kpu+ko-1
      il2(i)=npo(i-kpu+1)
      enddo
      do i=kpv,kpv+ko-1
      il2(i)=npo(i-kpv+1)
      enddo
      do i=kpo,kpo+ko-1
      il2(i)=npo(i-kpo+1)
      enddo
      do i=kpz,kpz+ko-1
      il2(i)=npo(i-kpz+1)
      enddo
      do i=kpt,kpt+ko-1
      il2(i)=npo(i-kpt+1)
      enddo
      do i=kpa,kpa+ko-1
      il2(i)=npo(i-kpa+1)
      enddo
      do i=kpr,kpr+ko-1
      il2(i)=npo(i-kpr+1)
      enddo
      do i=kpq,kpq+ko-1
      il2(i)=npo(i-kpq+1)
      enddo
      do i=kpo3,kpo3+ko-1
      il2(i)=npo(i-kpo3+1)
      enddo
      do i=kpcl,kpcl+ko-1
      il2(i)=npo(i-kpcl+1)
      enddo
      do i=kppr,kppr+ko-1
      il2(i)=npo(i-kppr+1)
      enddo
      do i=kpsun,kpsun+nsun-1
c     il2(i)=il2sun(i-kpsun+1)
      il2(i)=0
      enddo
c
      ienst=0
      iensi=0
      iens(1)=1
      iens(2)=ienst
      iens(3)=iensi
      iens(4)=1
      iens(5)=255
c
      call mpfdef(2,mpf)
c
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(lat1,lat2,iob2,io22,levs,si,sl,clat,fxs,ko,pokpa,fxp)
cmic$1 shared(jcap,nc,nctop,trig,ifax,eps,epstop,ss,sstop,slat)
cmic$1 shared(johf,nfldp,rf,mf)
cmic$1 shared(ksz,ksd,kst,ksq,kspsx,kspsy,ksu,ksv,ksps,kszs,kszsx,kszsy)
cmic$1 shared(kso3,kscl,kspr)
cmic$1 shared(kpu,kpv,kpo,kpz,kpt,kpr,kpsun,kslp)
cmic$1 shared(kpq,kpo3,kpcl,kppr)
cmic$1 shared(fxy,ntrace,ncloud)
cmic$1 shared(lpcl,lppr,lpo3)
cmic$1 private(lat,lan,rxs,oxs,zxs,zxi,tvx,n,i)
#endif
c
#ifdef SCM
      do lat=1,johf
#else
      do lat=2,johf
#endif
        do n=1,nfldp
          do i=1,io2
            fxp(i,n)=0.
          enddo
        enddo
#ifdef SCM
        call trss(nflds,ntrace,ncloud,ss,fxs)
#else
        call trss(nflds,ntrace,ncloud,trig,ifax,eps,epstop,ss,sstop,
     &            clat(lat),slat(lat),fxs(1,1    ))
#endif
#ifdef DBG
        do k=1,levs
          print 1960, k,fxs(1,kszs), fxs(1,ksps),
     +                fxs(1,ksd+k-1), fxs(1,ksz+k-1), fxs(1,ksu+k-1),
     +                fxs(1,kst+k-1), fxs(1,ksq+k-1)
1960      format(i5,7e13.5)
        enddo
#endif
        
      rmaxsiin=si(1)
      do k=1,levs+1
       rmaxsi=max(rmaxsi,si(k))
      enddo

      if (rmaxsi.gt.1. .or. rmaxsi.eq.0.) then
c
       rk  = rd_/cp_
       rk1 = rk + 1.
       rkinv=1./rk

       do i=1,io2
         do k=1,levs
           sihyb(i,k)=si(k)/fxs(i,ksps)/1000.+sl(k)
         enddo
         sihyb(i,levs+1)=0.
         do k=1,levs
           dif = sihyb(i,k)**rk1 - sihyb(i,k+1)**rk1
           dif = dif / (rk1*(sihyb(i,k)-sihyb(i,k+1)))
           slhyb(i,k) = dif**rkinv
c           slhyb(i,k)=0.5*(sihyb(i,k)+sihyb(i,k+1))
         enddo
       enddo
      else
        do i=1,io2
          do k=1,levs+1
            sihyb(i,k)=si(k)
          enddo
          do k=1,levs
            slhyb(i,k)=sl(k)
          enddo
        enddo
      endif

        call getrh(io2,io22,levs,slhyb,
     &             fxs(1,ksps    ),fxs(1,ksq    ),fxs(1,kst    ),
     &             rxs)
        call omega(io2,io22,levs,sihyb,slhyb,
#ifdef SCM
     &             fxs(1,ksu    ),
#else
     &             fxs(1,ksps    ),fxs(1,kspsx    ),fxs(1,kspsy    ),
     &             fxs(1,ksd    ),fxs(1,ksu    ),fxs(1,ksv    ),
#endif
     &             oxs)
        call hydro(io2,io22,levs,sihyb,slhyb,
     &             fxs(1,kszs    ),fxs(1,kst    ),fxs(1,ksq    ),
     &             zxs,zxi)
        call sig2p(io2,io22,levs,sihyb,slhyb,fxs(1,ksps    ),
#ifdef SCM
     &             fxs(1,ksd    ),fxs(1,ksz    ),oxs,
#else
     &             fxs(1,ksu    ),fxs(1,ksv    ),oxs,
#endif
     &             zxs,zxi,fxs(1,kst    ),rxs,fxs(1,ksq    ),
     &             fxs(1,kscl),fxs(1,kspr),fxs(1,kso3),
     &             ko,pokpa,lpcl,lppr,lpo3,
     &             fxp(1,kpu    ),fxp(1,kpv    ),fxp(1,kpo    ),
     &             fxp(1,kpz    ),fxp(1,kpt    ),fxp(1,kpr    ),
     &             fxp(1,kpq),
     &             fxp(1,kpcl),fxp(1,kppr),fxp(1,kpo3))
#ifdef DBG
        do k = 1, levs
            print 1964, k,fxp(1,kpu+k-1),fxp(1,kpv+k-1),fxp(1,kpo+k-1),
     +                  fxp(1,kpz+k-1), fxp(1,kpt+k-1), fxp(1,kpr+k-1),
     +                  fxp(1,kpq+k-1), fxs(1,ksps)
1964        format('PGB ',i5,8e13.5)
        enddo
#endif
        call sundry(io2,io22,levs,nsun,kslp,sihyb,
     &              fxs(1,kszs    ),fxs(1,ksps    ),
     &              fxs(1,ksq     ),fxp(1,kpz     ),
     &              fxp(1,kpsun    ))
        do n=1,nfldp
          do i=1,io2
            fxy(i,lat,n)=fxp(i,n)
          enddo
        enddo
      enddo
c
#ifndef SCM
c     call nntprt(fxy(1,1,1),io,jo,1./100.)
c
c  fill pole points with extrapolation
c
      do k=1,nfldp
        call polext(mpf(ipu(k)),io,fxy(1,2,k),fxy(1+io,2,k),
     &              fxy(1,1,k),fxy(1+io,1,k))
      call maxmin(fxy(1,1,k),io,jo,1,1,'after polext')
      enddo
c
c  smooth fields horizontally and compute absolute vorticity
c
      call llspc(io2,johf,ko,johf-2,ncpus,clat,slat,wlat,trig,ifax,
     &           fxy(1,1,kpu),fxy(1,1,kpv),fxy(1,1,kpo),
     &           fxy(1,1,kpz),fxy(1,1,kpt),fxy(1,1,kpa))
c
#endif
c
c  loop over groups of horizontal fields
c
      do k1=1,nfldp,ncpus
        k2=min(k1+ncpus-1,nfldp)
c
c  and round to the number of bits and engrib the field in parallel
c
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(k1,k2,pot,io2,johf,io,jo,mpf,nfldp,fxy,rf,mf)
cmic$1 shared(ipu,itl,il1,il2,icen,icen2,igen)
cmic$1 shared(idate,lenpds,jfhour,ids,mxbit,iens)
cmic$1 shared(grib,lgrib)
cmic$1 shared(rn)
cmic$1 private(k,kan,ierr,pok,poktop,lbm,ij)
#endif
        do k=k1,k2
          kan=k-k1+1
          lgrib(kan)=0
          print *,'k=',k,' ipu=',ipu(k)
          if(ipu(k).gt.0) then
            call rowsep(fxy(1,1,k))
cyh            call maxmin(fxy(1,1,k),io,jo,1,1,'fxy')
c
c         if(k.eq.1) then
c           call nntprt(fxy(1,1,k),io,jo,1./100.)
c         endif
#ifdef CRAY_THREAD
cmic$ guard
#endif
            call gribit(fxy(1,1,k),lbm,0,io,jo,mxbit,90.,
     &                  lenpds,2,icen,igen,0,
     &                  ipu(k),itl(k),il1(k),il2(k),
     &                  idate(4),idate(2),idate(3),idate(1),
     &                  1,jfhour,0,10,0,0,icen2,ids(ipu(k)),iens,
     &                  0.,0.,0.,0.,0.,0.,0.,0.,
     &                  grib(1,kan),lgrib(kan),ierr)
c
#ifdef CRAY_THREAD
cmic$ endguard
#endif
          endif
        enddo
c
c  write out grib messages sequentially
c
        do k=k1,k2
          kan=k-k1+1
          if((.not.lpo3).and.k.ge.kpo3.and.k.lt.kpcl) lgrib(kan) = 0
          if((.not.lpcl).and.k.ge.kpcl.and.k.lt.kppr) lgrib(kan) = 0
          if((.not.lppr).and.k.ge.kppr.and.k.lt.kpsun) lgrib(kan) = 0
          if(lgrib(kan).gt.0) then
            is=8
            lpds=ichar(grib(is+1,kan))*65536
     &          +ichar(grib(is+2,kan))*256
     &          +ichar(grib(is+3,kan))
            is=8+lpds
            lgds=ichar(grib(is+1,kan))*65536
     &          +ichar(grib(is+2,kan))*256
     &          +ichar(grib(is+3,kan))
            call wryte(npgb,lgrib(kan),grib(1,kan))
            print *,' grib1 written to ',npgb,' of length ',lgrib(kan)
          endif
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
