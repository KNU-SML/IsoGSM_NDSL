      subroutine surfeb(dtime,  month,
     2                   frTF,   pgcm,   Tair,   longw,
     3                 netsht,  vpair,    vpd,  rhoair,
     4                   wind,     ra,   binf,     dph,
     5                    qrt,    bkd,    sld,     wcr,
     6                    wpw,    smr,  dphnd,   smxnd,
     7                  expnd,  bubnd,  alpnd,   betnd,
     8                  gamnd,  vtype,    frt,   fflai,
     9                 sncvfr,   fswq,   frsn,    ftsf,
     a                   fcwt,   fcsn,   fsmc,    fsic,
     b                   ftnd,   latH,   senH,   grndH)

#include <define.h>
c
#ifdef VICLSM1
#include <vartyp.h>
#include <vic_surfeb.h>

************************************************************************
* source file:       surfeb.F
* PRGMMR:            Ji Chen   
* ORG:               ECPC/CRD/SIO/UCSD 
* DATE:              JUNE, JULY, AUGUST & SEPTEMBER 2003
* PRGM HISTORY:      modified from 4.0.3 UW VIC (calc_surf_energy_bal.c)
*
* ABSTRACT: calculates the surface temperature, in the
*           case of no snow cover.  Evaporation is computed using the
*           previous ground heat flux, and then used to comput latent 
*           heat in the energy balance routine.  Surface temperature
*           is found using the Frt Brent method (Numerical Recipies).
************************************************************************

* ------------------- input variables ----------------------------------
* -- model basic parameter
c
      integer month    ! current month
      real dtime       ! time step (second)

* -- atmosphere data

      real pgcm        ! pressure (pa)
      real Tair        ! air temperature (k)
      real rhoair      ! air density (kg/m3)
      real vpair       ! Actual vapor pressure of air (Pa)
      real vpd         ! vapor pressure deficit (Pa)
      real ra          ! aerodynamical resistance (s/m)
      real wind        ! wind speed (m/s)
      real netsht      ! net shortwave radiation (w/m2)
      real longw       ! longwave radiation (w/m2)

* -- soil parameters

      real binf        ! vic infiltration parameter (n/a)
      real dph(msl)    ! soil moisture layer thickness (m)
      real wcr(msl)    ! ~70% of field capacity (mm)
      real wpw(msl)    ! wilting point soil moisture (mm)
      real smr(msl)    ! residual moisture (mm)
      real sld(msl)    ! soil particle density (kg/m^3)
      real bkd(msl)    ! bulk density of soil layer (kg/m^3)
      real qrt(msl)    ! quartz content of soil (fraction)

      real dphnd(nsl)  ! soil thermal node thicknesses (m)
      real alpnd(nsl)  ! coef for computing soiltemp profile
      real betnd(nsl)  ! coef for computing soiltemp profile
      real gamnd(nsl)  ! coef for computing soiltemp profile
      real smxnd(nsl)  ! maximum soil moisture (m3/m3)
      real bubnd(nsl)  ! bubbling pressure of soil (cm)
      real expnd(nsl)  ! parameter for ksat with soil moisture (n/a)

      real kapnd(nsl)  ! soil layer thermal conductivity (w/m/k)
      real csnd(nsl)   ! soil layer heat capacity (j/m^3/k)
      real smnd(nsl)   ! soil moisture at node (m3/m3)
      real icend(nsl)  ! soil ice at node (m3/m3)

* -- vegetation parameters

      integer vtype    ! vegetation type
      real fflai       ! leaf area index (fraction)
      real frt(msl)    ! root content (fraction)

* -- land surface variables

      real sncvfr      ! snow cover fraction (fraction)
      real fcsn        ! snow on vegetation (m)
      real ftsf        ! snow surface temperature (k)
      real fswq        ! snow water equivalent (m)
      real frsn        ! snow density (kg/m3)
      real snowd       ! snow depth (m)

* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fcwt        ! dew and rain trapped on vegetation (m)
      real frtf        ! rain throughfall (m/time step)
      real fsmc(msl)   ! soil moisture (liquid + ice) (mm)
      real fsic(msl)   ! soil ice (mm)
      real ftnd(nsl)   ! soil tmperature profiles (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real lath        ! latent heat flux (w/m2)
      real senh        ! sensible heat flux (w/m2)
      real grndh       ! ground heat flux (w/m2)
*-----------------------------------------------------------------------

* ------------------- common block -------------------------------------
#include <vic_veglib.h>
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real d                ! vegetation displacement (m)
      real z0               ! surface roughness (m)
      real t0               ! land surface temperature (k)
      real d1               ! first soil layer thickness (m)
      real t_upper          ! up-boundary of temperature (k)
      real t_lower          ! low-boundary of temperature (k)

      real tsurf            ! surface temperature (k)
      real stabilityc       ! function name
      real root_brent       ! function name
      integer n, m          ! loop index

      real rhoh2o           ! water density (kg/m^3)
      real huge_resist
      real surf_dt
      parameter (surf_dt = 20.0,rhoh2o=1.e3, huge_resist = 1.e2)

      character*80  ctrfct  ! determine function
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,'vic debug -- enter surfeb',vtype,month
#endif

* ----------------------------------------------------------------------
* correct aerodynamic resistance for stability conditions
* ----------------------------------------------------------------------

c      if(vtype.gt.0) then
c         d  = veg_d(month, vtype)
c         z0 = veg_rough(month, vtype)
c      else
c         d  = 0
c         z0 = z0_soil
c      end if

      d = 0.0

c      if(fswq.gt.0.0) then
c         t0 = ftsf
c         z0 = 0.03
c      else

* because only ra is used for computing sensible and latent over
* no snow ground

         t0 = ftnd(1)
         z0 = 0.01

c      end if

#ifdef DBGVIC
      print *,'vic debug in surfeb before computing stability ',ra
#endif

      if (wind .gt. 0.0) then
         ra = ra/stabilityc(2.0, d, t0, tair, wind, z0)
      else
         ra = huge_resist
      end if  

#ifdef DBGVIC
      print *,'vic debug in surfeb after computing stability ',ra
#endif

* ----------------------------------------------------------------------
* find surface temperature using frt brent method
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,'vic debug in surfeb before node_soil'
      print *,'node_soil variables',msl,   nsl, dphnd,   dph,
     &             smxnd, expnd, bubnd,   sld,
     &               bkd,   qrt,  fsmc,  ftnd
#endif
      call node_soil(msl,   nsl, dphnd,   dph,
     &             smxnd, expnd, bubnd,   sld,
     &               bkd,   qrt,  fsmc,  ftnd,
     &              smnd, icend, kapnd,  csnd)

#ifdef DBGVIC
      print *, 'vic debug in surfeb after node_soil'
      print *,'node_soil variables',
     &              smnd, icend, kapnd,  csnd
#endif

* ----------------------------------------------------------------------
* calculate snow dph (h.b.h. 7.2.1)
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,'check frsn ',frsn
#endif

      if(sncvfr.gt.1.e-6) then
         snowd = rhoh2o * fswq / frsn / sncvfr ! in the units of m
      else
         snowd = 0.0
      endif

      d1 = dph(1)

* ----------------------------------------------------------------------
* added for temporary backwards compatability
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,'vic debug in surfeb before root_brent'
      print *,'dtime',dtime,pgcm,tair,rhoair,vpd,ra,wind,longw,netsht,
     1        binf,dph(1),wcr(1),wpw(1),smr(1),sld(1),bkd(1),
     2        qrt(1),frt(1),fsmc(1),fsic(1),dphnd(1),alpnd(1),
     3        betnd(1),gamnd(1),smxnd(1),bubnd(1),expnd(1),
     4        kapnd(1),csnd(1),smnd(1),icend(1),ftnd(1)
      print *,'vtpe ',vtype,fflai,frtf,fcwt,fcsn,ftsf,fswq,
     1        snowd
#endif
      t_upper = t0 + surf_dt
      t_lower = t0 - surf_dt

* ----------------------------------------------------------------------
* initialize vic_surfeb common block
* ----------------------------------------------------------------------

      edtime = dtime
      esncvfr= sncvfr
      epgcm  = pgcm
      etair  = tair
      erhoair= rhoair
      evpair = vpair
      evpd   = vpd
      era    = ra
      ewind  = wind
      elongw = longw
      enetsht= netsht
      ebinf  = binf
      
      do m = 1, msl
         edph(m) = dph(m)
         ewcr(m) = wcr(m)
         ewpw(m) = wpw(m)
         esmr(m) = smr(m)
         esld(m) = sld(m)
         ebkd(m) = bkd(m)
         eqrt(m) = qrt(m)
         efrt(m) = frt(m)
         efsmc(m) = fsmc(m)
         efsic(m) = fsic(m)
      end do
      do n = 1, nsl
         edphnd(n) = dphnd(n)
         ealpnd(n) = alpnd(n)
         ebetnd(n) = betnd(n)
         egamnd(n) = gamnd(n)
         esmxnd(n) = smxnd(n)
         ebubnd(n) = bubnd(n)
         eexpnd(n) = expnd(n)
         ekapnd(n) = kapnd(n)
         ecsnd(n)  = csnd(n)
         esmnd(n)  = smnd(n)
         eicend(n) = icend(n)
         eftnd(n) = ftnd(n)
      end do

      evtype = vtype
      efflai = fflai
      
      efrtf = frtf
      efcwt = fcwt
      efcsn = fcsn
      eftsf = ftsf
      efswq = fswq
      efrsn = frsn
      esnowd = snowd

      ctrfct = 'func_surfeb'
      tsurf = root_brent(t_upper, t_lower, ctrfct)

      fcwt = mfcwt
      frtf = mfrtf
      
      do m = 1, msl
         fsmc(m) = mfsmc(m)
      end do
      do n = 1, nsl
         ftnd(n) = mftnd(n)
      end do

      lath = olath
      senh = osenh
      grndh = ogrndh

#ifdef DBGVIC
      print *,'vic debug -- end of surfeb  tsurf (k) ',tsurf
#endif
#endif
      return
      end
