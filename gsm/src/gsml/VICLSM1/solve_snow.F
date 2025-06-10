      subroutine solve_snow(msl,     dtime,    month,     rfall,
     2                    sfall,      pgcm,     tair,     longw, 
     3                  netshts,   netshtv,   rhoair,     vpair,
     4                      vpd,        ra,        u,       wcr,
     5                      wpw,     vtype,      ost,       frt,
     6                    fflai,     tgrnd,     fsmc,      fsic,
     7                     fswq,      frsn,     ftsf,      ftpk,
     8                     fsfw,      fpkw,     fcwt,      fcsn,
     9                     frtf,    fsnowm,   sncvfr,     laths,
     a                    senhs)

#include <vartyp.h>

************************************************************************
* source file:       solve_snow.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (solve_snow.c)
*
* abstract: this routine is to handle the various calls for 
*           solving the various components of the vic snow scheme.
************************************************************************

* ------------------- input variables ----------------------------------

* -- model basic parameters

      integer msl           ! number of soil layer
      integer month         ! current month

      real dtime            ! time step (s)

* -- atmosphere variables

      real rfall            ! rainfall (m/time step)
      real sfall            ! snowfall (m/time step)
      real pgcm             ! pressure (pa)
      real tair             ! air temperature (k)
      real longw            ! long wave (w/m2)
      real netshts          ! net short wave radiation on snow ground (w/m2)
      real netshtv          ! net short wave radiation no snow ground (w/m2)
      real vpair            ! actual vapor pressure of air (pa)
      real vpd              ! vapor pressure deficit (pa)
      real rhoair           ! air density (kg/m3)
      real u(3)             ! wind speed (m/s)
      real ra(3)            ! aerodynamical resistance (s/m)

* -- soil parameters

      real wcr(msl)         ! ~70% of field capacity (mm)
      real wpw(msl)         ! wilting point soil moisture (mm)

* -- vegetation parameters

      integer vtype         ! vegetation class
      integer ost           ! vegetation overstory (1: with, 0: no)
      real frt(msl)         ! root content (fraction)
      real fflai            ! leaf area index

* -- land surface variable

      real tgrnd            ! ground surface temperature (k)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fswq             ! snow water equivalent (m h2o)
      real frsn             ! snow density (kg/m3)
      real ftsf             ! snow surface temperature (k)
      real ftpk             ! snow pack temperature (k)
      real fsfw             ! snow surface water equivalent (m h2o)
      real fpkw             ! snow pack water equivalent (m h2o)
      real fcwt             ! canopy intercepted water (m h2o)
      real fcsn             ! canopy intercepted snow (m h2o)
      real fsmc(msl)        ! soil moisture (liquid + ice) (mm)
      real fsic(msl)        ! soil ice (mm)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real frtf             ! rain through fall (m/time step)
      real fsnowm           ! snow melt (m/time step)
      real sncvfr           ! snow cover fraction (fraction)
      real laths            ! latent heat flux (w/m2)
      real senhs            ! sensible heat flux (w/m2)
* ----------------------------------------------------------------------

* ------------------- common block -------------------------------------
#include <vic_veglib.h>
* ----------------------------------------------------------------------

* ------------------- local variables ----------------------------------
      real fstf             ! snow through fall (m/time step)
      real cnplongw         ! longwave radiation to canopy (w/m2)
      real longout          ! land emitted longwave radiation (w/m2)
      real netsht           ! net short wave radiation (w/m2)

      real vpcnp            ! canopy air saturated vapor pressure
      real le, ls           ! latent heat of vaporization/sublimation (j/kg)
      real lathsc           ! latent heat fluxes over canopy (w/m2)
      real senhsc           ! sensible heat fluxes over canopy (w/m2)
      real evapc            ! canopy evapotranspiration (w/m2)
      real snowd            ! snowpack depth (m)
      real rad              ! land surface net radiation (w/m2)
      real sfat             ! radiation attentuation factor

      real old_fswq         ! fswq before computing snow melt (m)
      real tmptsf4          ! ftsf**4
      real tmptgr4          ! tgrnd**4
      real tmp1             ! sncvfr*sigma

      real mfcwt            ! temporary canopy intercepted water (m h2o)
      real mfsmc(msl)       ! temporary soil moisture (liquid + ice) (mm)
      integer m             ! loop index

      real canopy_evap      ! function name
      real snow_density     ! function name

      real mxcvswq          ! minimum swq for fully covering a grid (m)
      parameter (mxcvswq=0.00)

      real rhoh2o           ! water density (kg/m^3)
      real t0c              ! ice/water mix temperature (k)
      real sigma            ! stefan-boltzmann constant (w/m2/k-4)
      parameter (t0c=273.15, sigma=5.67e-8,rhoh2o=1.e3)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* initialize the output and other variables
* ----------------------------------------------------------------------

      lathsc = 0.0
      senhsc = 0.0
      evapc  = 0.0
      fstf   = 0.0

* ----------------------------------------------------------------------
* check for thin snowpack which only partially covers grid cell
* ----------------------------------------------------------------------

      if(fswq .lt. mxcvswq) then
         sncvfr = fswq/mxcvswq
      else 
         sncvfr = 1.0
      endif

* ----------------------------------------------------------------------
* if vegetation overstory exists
* ----------------------------------------------------------------------

      if(ost.eq.1) then

* ----------------------------------------------------------------------
* compute the snow interceptation over whole grid
* ----------------------------------------------------------------------

         tmp1 = sncvfr*sigma
         tmptsf4 = ftsf**4
         tmptgr4 = tgrnd**4

         cnplongw = longw + tmp1*tmptsf4 + (sigma-tmp1)*tmptgr4

         if(fcsn.gt.0.or.sfall.gt.0) then
            call snow_intcept(dtime,    fflai, rhoair,   tair,
     &                        vpair,      vpd,  ra(2),   u(2),
     &                      netshts, cnplongw,  rfall,  sfall,
     &                         pgcm,     fcwt,   fcsn,   frtf,
     &                         fstf,   lathsc, senhsc)

         else                            ! canopy no snow

* ----------------------------------------------------------------------
* calculate the net radiation at the canopy surface, using the canopy 
*    temperature. the outgoing longw is subtracted twice, because the 
*    canopy radiates in two directions
* ----------------------------------------------------------------------

            longout = sigma * tair**4
            rad   = netshtv + cnplongw - 2*longout

            evapc = canopy_evap(msl,  vtype,   dtime,   fflai,
     &                         fcwt,  rfall,    fsmc,    fsic,
     &                          wcr,    wpw,     frt,     rad,
     &                       rhoair,  vpair,     vpd,   ra(3),
     &                      netshtv,   tair,    pgcm,   mfcwt,
     &                         frtf,  mfsmc)

            fcwt = mfcwt
            do m = 1, msl
               fsmc(m) = mfsmc(m)
            enddo

            if(tair.gt.t0c) then
               le = (2.501 - 0.002361 * (tair-t0c)) * 1.0e6
               evapc = le*evapc*rhoh2o
            else
               ls = (677.-0.07*(tair-t0c))*4.1868*1000.0
               evapc = ls*evapc*rhoh2o
            end if

         end if
      else
         frtf = rfall
         fstf = sfall
      end if            

#ifdef DBGVIC
      print *,'in solve_snow lathsc/evapc=',lathsc,evapc
      print *,'in solve_snow       senhsc=',senhsc
#endif

      if((fswq+fstf).gt.0.0)then
         
* ----------------------------------------------------------------------
* check for thin snowpack which only partially covers grid cell
* the snow melt is modeled over the sncvfr
* ----------------------------------------------------------------------

         if((fswq+fstf) .lt. mxcvswq) then
            sncvfr = (fswq+fstf)/mxcvswq
         else
            sncvfr = 1.0
         endif

* ----------------------------------------------------------------------
* shortwave radiation attentuation factor because of canopy
* ----------------------------------------------------------------------

         if(ost.eq.1)then
#define CORR
#ifdef CORR
            sfat = exp(-0.1*veg_sfat(vtype)*fflai)
#else
            sfat = exp(-veg_sfat(vtype)*fflai)
#endif
         else
            sfat = 1.0
         endif
         netsht = netshts*sfat

* ----------------------------------------------------------------------
* compute snow melt over the snow_coverage_fraction (sncvfr)
* convert the snow variables to the fraction of sncvfr
* ----------------------------------------------------------------------

         fstf = fstf/sncvfr
         fswq = fswq/sncvfr
         fsfw = fsfw/sncvfr
         fpkw = fpkw/sncvfr

         old_fswq = fswq

* ----------------------------------------------------------------------
* calculate snow depth (h.b.h. 7.2.1) 
* ----------------------------------------------------------------------

         snowd = 1000. * fswq / frsn        ! in the units of m

* ----------------------------------------------------------------------
* call snow pack accumulation and ablation algorithm
* ----------------------------------------------------------------------
         call snow_melt(dtime,   vtype,   month,  rhoair,
     &                  ra(3),    u(3),   vpair,     vpd,
     &                   tair,   tgrnd,  netsht,   longw,
     &                   pgcm,    frtf,    fstf,   snowd,
     &                   frsn,    fswq,    ftsf,    ftpk,
     &                   fsfw,    fpkw,  fsnowm,   laths,
     &                  senhs)

         laths = laths*sncvfr + lathsc + evapc
         senhs = senhs*sncvfr + senhsc

#ifdef DBGVIC
         print *,'in solve_snow laths=',laths,lathsc,evapc
         print *,'in solve_snow senhs=',senhs,senhsc
#endif

* ----------------------------------------------------------------------
* compute snow parameters 
* ----------------------------------------------------------------------

         if(fswq .gt. 0.) then

#ifdef DBGVIC
            print *,'in solve_snow befsndens',frsn,fstf,old_fswq,snowd
#endif

* ----------------------------------------------------------------------
* calculate snow density
* ----------------------------------------------------------------------

            frsn = snow_density(fstf, tair, old_fswq, snowd,
     &                         dtime, ftsf)

#ifdef DBGVIC
            print *,'in solve_snow after snow_density',frsn
#endif
         end if

* ----------------------------------------------------------------------
* compute snow melt over the snow_coverage_fraction (sncvfr)
* convert the snow parameters to sncvfr
* ----------------------------------------------------------------------

         fswq  = fswq*sncvfr
         fsfw  = fsfw*sncvfr
         fpkw  = fpkw*sncvfr
         fsnowm= fsnowm*sncvfr

      else

* ----------------------------------------------------------------------
* ground snow not present
* ----------------------------------------------------------------------

         laths = lathsc + evapc
         senhs = senhsc
      end if

      if(fswq.le.1.e-6) then
         frsn = 50.0
         fswq =  0.0
         fsfw =  0.0
         fpkw =  0.0
         ftsf = min(tair,t0c)
         ftpk = min(tgrnd,t0c)
      end if

      return
      end
