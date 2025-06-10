#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
#define JCAP jcap
#define LATG latg
#define NPES npes
#include <define.h>
      subroutine lwrchou(ipts,hlw,slwnc,tlwuc,
#ifdef CLR
     1              slwn0,tlwu0,
#endif
     1              pl2,temp,rh2o,qo3,taulw,cld0)
c
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
#include <rdparm.h>
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#else
#include <commpi.h>
#endif
#endif
c ---------------------------------------------------------------------
c [06/07/00] fanglin yang
c     1) add "save" statement for data blocks.
c     2) change all "real*8" in longwave subroutine to "real".
c     3) fix nt=1 in irrad1.F, and rename h11,h21,h81,c1 and o1
c        as hh11,hh21,hh81,cc1 and oo1 to avoid memeory confusion
c        for shared-memeroy computers
c ---------------------------------------------------------------------
c [05/01/00] fanglin yang
c     computes ir fluxes due to water vapor, co2, o3, trace gases
c     (n2o, ch4, cfc11, cfc12, cfc22, co2-minor),clouds, and aerosols.
c --input variables
c     pl2    : sigma-layer pressure (mb), k=1-->toa
c     temp   : integer-layer temperature (k),k=1--> top layer,
c              with temp(lp1) = surface temperature
c     rh2o   : water vapor mixing ratio (g/g), k=1-->top layer
c     qo3    : o3 mixing ratio (g/g), k=1-->top layer
c     taulw  : cloud optical depth, k=1-->top layer
c     cld0   : cloud fraction, k=1-->top layer
c --output variables
c     hlw    : longwave cooling rate(k/day), negative
c     slwnc  : surface net lw flux (upward +) (1.e-3 w/m^2)
c     tlwuc  : top upward(+) lw flux (1.e-3 w/m^2)
c     slwn0  : clear-sky surface net lw flux(upward +) (1.e-3 w/m^2)
c     tlwu0  : clear-sky top upward(+) lw flux (1.e-3 w/m^2)
c ---------------------------------------------------------------------
c
c -- input
      dimension pl2(imax,lp1), temp(imax,lp1), rh2o(imax,l),
     >       qo3(imax,l),   taulw(imax,l), cld0(imax,l)
c -- output, radiative fluxex are in cm.g.s. system
      dimension hlw(imax,l), slwnc(imax), tlwuc(imax),
     >       slwn0(imax), tlwu0(imax)
c
c --temporary, all radiative fluxes are in w/m^2, downward(+)
      dimension pir(imax,lp1), ta(imax,l), wa(imax,l), oa(imax,l),
     >       ccld(imax,l),  taucl(imax,l), taual(imax,l,10),
     >       tg(imax)
      dimension flxd(imax,lp1), flxu(imax,lp1), flx(imax,lp1),
     >       flxdc(imax,lp1),flxuc(imax,lp1),flc(imax,lp1),
     >       dfdts(imax,lp1),st4(imax)
      real n2oppm, n2o
      save first
      logical first
      data first/.true./
c
c * if model top is above 30mb, please use high=.true.
      logical high,trace
cc    data high/.true./
cc    data trace/.false./
       high=.true.
       trace=.false.  !without ch4, cfcs, and n2o

c   level index separating high and middle clouds (ict)
c   level index separating middle and low clouds   (icb)
        ict       =15    !about 440mb for mls standard atmos.
        icb       =19    !about 702mb for mls standard atmos.
#ifndef CO2PPM
        co2ppm    =348.0
#else
        co2ppm    =CO2PPM
#endif
#ifndef CH4PPM
        ch4ppm    =1.75
#else
        ch4ppm    =CH4PPM
#endif
#ifndef N2OPPM
        n2oppm    =0.28
#else
        n2oppm    =N2OPPM
#endif
        f22ppm    =0.2e-03      !ppmv
        f11ppm    =0.2145E-3    !cfc11 concentration in ppm
        f12ppm    =0.3711E-3    !cfc12 concentration in ppm
c -----exchange unit
        co2       =co2ppm*1.e-6
        ch4       =ch4ppm*1.e-6
        n2o       =n2oppm*1.e-6
        cfc11     =f11ppm*1.e-6
        cfc12     =f12ppm*1.e-6
        cfc22     =f22ppm*1.e-6
c
c -- exchange temporary variables
      do i=1,ipts
      do k=1,l
         pir(i,k)     =pl2(i,k)
         ta(i,k)      =temp(i,k)
         wa(i,k)      =rh2o(i,k)
         oa(i,k)      =qo3(i,k)
         ccld(i,k)    =cld0(i,k)
         taucl(i,k)   =taulw(i,k)
         do kk=1,10
         taual(i,k,kk)=0.    !aerosol optical depth, set to zero
         enddo
      enddo
         tg(i)        =temp(i,lp1)
         pir(i,lp1)  =pl2(i,lp1)
      enddo
      if (first) then
#ifndef NOPRINT
#ifdef MP
      if( mype.eq.master ) then
#endif
        print *,'lwrchou: co2 concentration is ',co2
        print *,'lwrchou: ch4 concentration is ',ch4
        print *,'lwrchou: n2o concentration is ',n2o
#ifdef MP
      endif
#endif
#endif
        first=.false.
      endif
         

c
c     ------------------------------------------------------------
c     ------------------------------------------------------------
      call irrad1 (ipts,imax,l,pir,ta,wa,oa,tg,co2,
     *     n2o,ch4,cfc11,cfc12,cfc22,
     *     taucl,ccld,ict,icb,taual,
     *     high,trace,flxd,flxu,flx,flxdc,flxuc,flc,dfdts,st4)
c     ------------------------------------------------------------
c     ------------------------------------------------------------
c     constant for calculating heating rate (k/day)
c     hrate=daysecond*gravity/cpc=24*3600*9.80616/1004.0
c     where cpc is the specific heat of dry air (j/deg/kg)
      hrate = 843.8767171

      do i=1,ipts
      do k=1,l
       hlw(i,k)=0.01*hrate*(flx(i,k)-flx(i,k+1))/(pir(i,k+1)-pir(i,k))
      enddo
       slwnc(i) =-1000.*flx(i,lp1)
       slwn0(i) =-1000.*flc(i,lp1)
       tlwuc(i) =-1000.*flxu(i,1)
       tlwu0(i) =-1000.*flxuc(i,1)
      enddo

      return
      end
