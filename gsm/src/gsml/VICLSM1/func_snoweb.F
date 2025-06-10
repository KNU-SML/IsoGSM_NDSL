      function func_snoweb(tsf)

#include <vartyp.h>

************************************************************************
* source file:       func_snoweb.f
* purpose:           calculate snow pack energy balance
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (snowpackenergybalance)
*
* notes:
*    in c language, the values of function parameters are passed one-way, 
*    i.e., the value is from the upper to down; but two-way in fortran
*
* reference: bras, r. a., hydrology, an introduction to hydrologic
*            science, addisson wesley, inc., reading, etc., 1990.
************************************************************************

* ------------------- input variables ----------------------------------
      real tsf            ! snow surface temperature (c)
* ----------------------------------------------------------------------

* ------------------- output variables ----------------------------------
      real func_snoweb    ! residual snow energy (w/m2)
* ----------------------------------------------------------------------

* ------------------- common blocks ------------------------------------
#include <vic_snoweb.h>
#include <vic_veglib.h>
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real le              ! latent heat of vaporization (j/kg)
      real ls              ! latent heat of sublimation (j/kg)
      real d               ! vegetation displace (m)
      real z0              ! vegetation roughness (m)
      real rac             ! corrected resistance (s/m)
      real longradout      ! emitted longwave radiation (w/m2)
      real rad             ! net radiation (w/m2)
      real vpsnow          ! saturated vapor pressure in snow pack (pa)
      real vapormf         ! mass flux of vapor to or from snow (m/s)
      real advectede       ! energy advected by precipitation (w/m2)
      real deltcc          ! change in cold content (w/m2)
      real spkh            ! snow pack heat flux (w/m2)
      real restterm        ! rest term in surface energy balance (w/m2)
      real stabilityc, svp ! function names

      real cp              ! specific heat of moistair at constant (j/kg/k)
      real lf              ! latent heat of freezing (j/kg) at 0c
      real rhoh2o          ! water density (kg/m^3)
      parameter (cp=1013.0,lf=3.337e5,rhoh2o=1.e3)

      real sigma           ! stefan-boltzmann constant (w/m^2/k^-4)
      real t0c             ! ice/water mix temperature (k)
      real huge_resist
      parameter (huge_resist=1.e2, sigma=5.67e-8, t0c=273.15)

      real ch_water        ! volumetric heat capacity (j/(m3*c)) of water
      real ch_ice          ! volumetric heat capacity (j/(m3*c)) of ice
      parameter (ch_water=4186.8e3,ch_ice=2100.0e3)
* ----------------------------------------------------------------------

      le = (2.501 - 0.002361 * tsf) * 1.0e6
      ls = (677. - 0.07 * tsf) * 4.1868 * 1000.0
c      d  = veg_d(emonth, evtype)
c      z0 = veg_rough(emonth, evtype)

* ----------------------------------------------------------------------
* correct aerodynamic conductance for stable conditions
* ----------------------------------------------------------------------
      d = 0.0
      z0 = 0.03

      if (ewind .gt. 0.0) then
         rac = era/stabilityc(2.0, d, tsf+t0c, etair+t0c, ewind, z0)
      else
         rac = huge_resist
      end if

* ----------------------------------------------------------------------
* calculate longwave exchange and net radiation
* ----------------------------------------------------------------------
      longradout = sigma * (tsf+t0c)**4

      rad = enetsht + elongw - longradout
  
* ----------------------------------------------------------------------
* calculate the sensible heat flux
* ----------------------------------------------------------------------
      osenhs = erhoair * cp * (etair - tsf)/rac

#ifdef DBGVIC
      print *,'in func_snoweb osenhs',erhoair,cp,etair,tsf,era,rac
#endif
* ----------------------------------------------------------------------
* calculate the mass flux of ice to or from the surface layer
* calculate saturated vapor pressure in snow pack, (eq. 3.32, bras 90)
* ----------------------------------------------------------------------
      vpsnow  = svp(tsf+t0c)
      vapormf = erhoair * (0.622/epgcm) * (evpair - vpsnow)/rac
      vapormf = vapormf/rhoh2o

#ifdef DBGVIC
      print*,'in func_snoweb vapormf ',vapormf,erhoair,epgcm,evpair,
     1       vpsnow,rac,enetsht,elongw,longradout
#endif

      if (evpd .eq. 0.0 .and. vapormf .lt. 0.0) vapormf=0.0
  
* ----------------------------------------------------------------------
* calculate latent heat flux
* ----------------------------------------------------------------------

      if (tsf .ge. 0.0) then

* ----------------------------------------------------------------------
* melt conditions: use latent heat of vaporization 
* ----------------------------------------------------------------------
         olaths = le * vapormf * rhoh2o
      else 

* ----------------------------------------------------------------------
* accumulation: use latent heat of sublimation (eq. 3.19, bras 1990)
* ----------------------------------------------------------------------
         olaths = ls * vapormf * rhoh2o
      end if
  
* ----------------------------------------------------------------------
* calculate advected heat flux from rain
* ----------------------------------------------------------------------
      advectede = (ch_water * (etair-tsf) * efrtf) / edtime
  
* ----------------------------------------------------------------------
* calculate change in cold content
* ----------------------------------------------------------------------
      deltcc = ch_ice * esfswq * (eoldtsf-tsf)/ edtime

* ----------------------------------------------------------------------
* discussion: (ji 2/2004) the contribution from deltcc to the snow 
*    surface energy balance should be not very high and we should limit
*    the valud of deltcc. the reason is that in the energy balance equation
*    all the energy terms are facing to a virtual snow surface instead
*    of the surface snow layer.
* ----------------------------------------------------------------------

c      deltcc = deltcc/10.0
c      if(abs(deltcc).gt.100) deltcc = sign(100.0, deltcc)
c      deltcc = 0.0

#ifdef DBGVIC
      print *,'in func_snoweb ---',ch_ice,esfswq,eoldtsf,tsf,edtime
#endif

* ----------------------------------------------------------------------
* calculate snow pack heat flux from the ground surface
* ----------------------------------------------------------------------
      if(esnowd.gt.0.) then
         spkh = (2.9302e-6)* efrsn*efrsn * (etgrnd-tsf)/esnowd
      else 
         spkh = 0.0
      end if

c      if(abs(spkh).gt.100) spkh = sign(100.0, spkh)

* ----------------------------------------------------------------------
* calculate net energy exchange at the snow surface
* ----------------------------------------------------------------------
#ifdef DBGVIC
      print *,'in func_snoweb +++',tsf,rad,osenhs,olaths,advectede,
     1        spkh, deltcc, edtime, esfwater,lf,etgrnd,esnowd
#endif

      restterm  = rad + osenhs + olaths + advectede + spkh + deltcc
      orefreeze = (esfwater*lf*rhoh2o)/edtime

#ifdef DBGVIC
      print *,'in func_snoweb ===',restterm, orefreeze
#endif

* ----------------------------------------------------------------------
* please check wigmosta, vail and lettenmaier (p 1670, wrr 1994)
*   for the computation of snow melt energy
*   bascially, when tsf < 0c no snow melt, 
*                 tsf = 0, allowing the heat flux avaiable for snowmelt
* ----------------------------------------------------------------------

      if (abs(tsf).lt.1.e-6 .and. restterm .ge. (-orefreeze)) then
         orefreeze = -restterm           ! available energy input over
                                         ! cold content used to melt
         restterm  = 0.0
      else
         restterm  = restterm + orefreeze  ! add this value to the pack
      end if

      func_snoweb  = restterm

#ifdef DBGVIC
      print *,'in func_snoweb ***',restterm, orefreeze
#endif

      return
      end
