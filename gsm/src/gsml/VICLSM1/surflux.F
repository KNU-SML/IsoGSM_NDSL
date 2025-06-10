      subroutine surflux(msl,    nsl,  dtime,  month,
     2                 rfall,  sfall,   pgcm,   tgcm,
     3                 flwds,   sols,  vpair,    vpd,
     4                rhoair,      u,     ra,   binf,
     5                   dph,    qrt,    bkd,    sld,
     6                   wcr,    wpw,    smr,  dphnd,
     7                 smxnd,  expnd,  bubnd,  alpnd,
     8                 betnd,  gamnd,  vtype,    frt,
     9                 fflai,  lstsn,   fswq,   frsn,
     a                  ftsf,   ftpk,   fsfw,   fpkw,
     b                  fcwt,   fcsn,   fsmc,   fsic,
     c                  ftnd,  falbd,   lath,   senh,
     d                 grndh,   frtf, fsnowm, fsnowe,
     e                sncvfr)

#include <vartyp.h>

************************************************************************
* source file:       surface_fluxes.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (surface_fluxes.c)
*
* abstract: this routine computes surface fluxes, and solves 
*           the snow accumulation and ablation algorithm. solutions
*           are for the current snow band and vegetation type
************************************************************************

* ------------------- input variables ----------------------------------

* -- model basic parameters

      integer msl           ! number of soil layer
      integer nsl           ! number of soil thermal nodes
      integer month         ! current month

      real dtime            ! time step (second)

* -- atmosphere variables

      real rfall            ! rainfall (m/time step)
      real sfall            ! snowfall (m/time step)
      real pgcm             ! pressure (pa) (vic: use elev(m) to get p)
      real tgcm             ! air temperature (k)
      real flwds            ! long wave (w/m^2)
      real sols             ! short wave (w/m^2)
      real vpair            ! actual vapor pressure of air (pa)
      real vpd              ! vapor pressure deficit (pa)
      real rhoair           ! air density (kg/m^3)
      real u(3)             ! wind speed (m/s)
      real ra(3)            ! aerodynamical resistance (s/m)

* -- soil parameters

      real binf             ! vic infiltration parameter
      real dph(msl)         ! soil layer thickness (m)
      real qrt(msl)         ! quartz content of soil (fraction)
      real bkd(msl)         ! bulk density of soil layer (kg/m^3)
      real sld(msl)         ! soil particle density (kg/m^3)
      real wcr(msl)         ! ~70% of field capacity (mm)
      real wpw(msl)         ! wilting point soil moisture (mm)
      real smr(msl)         ! residual moisture (mm)

      real dphnd(nsl)       ! soil node thickness (m)
      real smxnd(nsl)       ! maximum soil moisture at soil node (m3/m3)
      real expnd(nsl)       ! parameter for ksat with soil moisture (n/a)
      real bubnd(nsl)       ! bubbling pressure of soil (cm)
      real alpnd(nsl)       ! coef for computing soiltemp profile
      real betnd(nsl)       ! coef for computing soiltemp profile
      real gamnd(nsl)       ! coef for computing soiltemp profile

* -- vegetation parameters

      integer vtype         ! vegetation class
      real frt(msl)         ! root content (fraction)
      real fflai            ! leaf area index

* -- land surface variable

      integer lstsn         ! surface snow age (number of model time step)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fswq             ! snow water equivalent (m h2o)
      real frsn             ! snow density (kg/m^3)
      real ftsf             ! snow surface temperature (k)
      real ftpk             ! snow pack temperature (k)
      real fsfw             ! snow surface water equivalent (m h2o)
      real fpkw             ! snow pack water equivalent (m h2o)
      real fcwt             ! canopy intercepted water (m h2o)
      real fcsn             ! canopy intercepted snow (m h2o)
      real fsmc(msl)        ! soil moisture (liquid + ice) (mm)
      real fsic(msl)        ! soil ice (mm)
      real ftnd(nsl)        ! soil temperature profile (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real falbd            ! land surface albedo (fraction)
      real lath             ! latent heat flux (w/m^2) (+: to surface)
      real senh             ! sensible heat flux (w/m^2) (+: to surface)
      real grndh            ! ground heat flux (w/m^2) (+: to surface)
      real frtf             ! rain throughfall (m/time step)
      real fsnowm           ! snow melt (m/time step)
      real fsnowe           ! snow evap (w/m2)
      real sncvfr           ! snow cover fraction (fraction)
* ----------------------------------------------------------------------

#include <vic_veglib.h>

* ---------------------- local variables -------------------------------
      real falbds           ! snow related surface albedo (fraction)
      real falbdv           ! without snow surface albedo (fraction)
      real tgrnd            ! ground surface temperature (k)
      real laths            ! latent heat flux over snow surface (w/m^2)
      real senhs            ! sensible heat over snow surface (w/m^2)
      real netshts          ! net shortwave radiation for snow ground (w/m2)
      real netshtv          ! net shortwave radiation without snow (w/m2)
      integer ost           ! vegetation overstory (1: with, 0: no)

      real snowage          ! snow age (day)
      integer secpday
      parameter (secpday=86400)
      real snow_albedo      ! function name
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,'vic debug -- enter surflux',fsmc,fsic
#endif

* ----------------------------------------------------------------------
* initialize variables
* ----------------------------------------------------------------------

      ost = veg_ost(vtype)
      lath  = 0.0
      senh  = 0.0
      grndh = 0.0
      frtf  = 0.0
      fsnowm = 0.
      fsnowe = 0.

      tgrnd  = ftnd(1)       ! ground surface temperature
      laths  = 0.
      senhs  = 0.
      sncvfr = 0.

      falbds = 0.85          ! new snow albedo

      falbdv = veg_alb(month, vtype)

      if(fswq.gt.0.or.sfall.gt.0.0)then
* ----------------------------------------------------------------------
* compute snow pack albedo 
* ----------------------------------------------------------------------
         if(lstsn.lt.0) then            ! due to snow existing on ground
            print *,'warning in surflux for setting lstsn =',lstsn
            print *,' in surflux fswq, sfall',fswq,sfall
            lstsn = 0
         endif
         
         snowage = float(lstsn)*dtime/float(secpday)
         falbds  = snow_albedo(snowage, ftsf)
      end if

      netshts = (1. - falbds)*sols
      netshtv = (1. - falbdv)*sols

* ----------------------------------------------------------------------
* solve for overstory canopy hydrology, snow pack accumulation & ablation
* ----------------------------------------------------------------------

      if(fswq.gt.0.or.sfall.gt.0.or.(fcsn.gt.0.and.ost.eq.1)) then

#ifdef DBGVIC
         print *,'vic debug in surflux before solve_snow'
         print *,' fswq=',fswq,' fcsn=',fcsn,' sfall=',sfall
         print *,' bef solve_snow fsmc,fsic',fsmc,fsic
#endif

         call solve_snow(msl,     dtime,    month,     rfall,
     2                 sfall,      pgcm,     tgcm,     flwds, 
     3               netshts,   netshtv,   rhoair,     vpair,
     4                   vpd,        ra,        u,       wcr,
     5                   wpw,     vtype,      ost,       frt,
     6                 fflai,     tgrnd,     fsmc,      fsic,
     7                  fswq,      frsn,     ftsf,      ftpk,
     8                  fsfw,      fpkw,     fcwt,      fcsn,
     9                  frtf,    fsnowm,   sncvfr,     laths,
     a                 senhs)

         fsnowe = - laths   ! positive upward for fsnowe (w/m2)

#ifdef DBGVIC
         print *,' after solve_snow fsmc,fsic',fsmc,fsic
#endif
      else

#ifdef DBGVIC
         print *,'vic in surflux no snow cover '
#endif
         frtf = rfall

      end if

      falbd = falbds*sncvfr + (1.-sncvfr)*falbdv

* ----------------------------------------------------------------------
* solve energy balance components for ground
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *, 'vic debug in surflux before surfeb ra=',ra
      print *,'in surflux fflai',fflai,' vtype',vtype
#endif

      call surfeb(dtime,  month,
     2             frtf,   pgcm,   tgcm,  flwds,
     3          netshtv,  vpair,    vpd, rhoair,
     4             u(1),  ra(1),   binf,    dph,
     5              qrt,    bkd,    sld,    wcr,
     6              wpw,    smr,  dphnd,  smxnd,
     7            expnd,  bubnd,  alpnd,  betnd,
     8            gamnd,  vtype,    frt,  fflai,
     9           sncvfr,   fswq,   frsn,   ftsf,
     a             fcwt,   fcsn,   fsmc,   fsic,
     b             ftnd,   lath,   senh,  grndh)

#ifdef DBGVIC
      print *, 'in surflux lath=',lath,' laths=',laths
      print *, 'in surflux senh=',senh,' senhs=',senhs
#endif

      lath = (1.0-sncvfr)*lath + laths
      senh = (1.0-sncvfr)*senh + senhs

      if(abs(lath).gt.800) then
         print *,'** warning for lath in surflux ',lath, laths, sncvfr
      endif
      if(abs(senh).gt.800) then
         print *,'** warning for senh in surflux ',senh, senhs, sncvfr
c         stop
      endif

#ifdef DBGVIC
      print *, 'vic debug -- end of surflux'
#endif

      return
      end
