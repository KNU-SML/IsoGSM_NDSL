      function func_surfeb(ts)

#include <define.h>
c
#ifdef VICLSM1
#include <vartyp.h>
#include <vic_surfeb.h>

************************************************************************
* source file:       func_surfeb.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (func_surf_energy_bal.c)
*
* abstract: this subroutine computes the surface energy balance for bare 
*   soil and vegetation uncovered by snow.  it computes outgoing longwave,
*   sensible heat flux, ground heat flux, and storage of heat in the thin
*   upper layer, based on the given surface temperature.
*        the energy balance equation used comes from xu liangs paper 
*  "insights of the ground heat flux in land surface 
*   parameterization schemes."
************************************************************************

* ------------------- input variables ----------------------------------
      real ts              ! surface temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ----------------------------------
      real func_surfeb     ! surface energy balance error (w/m2)
* ----------------------------------------------------------------------

* ------------------- common blocks ------------------------------------
#include <vic_veglib.h>
      integer kdt,jdt,inistp,limlow,maxstp,numsum,nummax
      common/comfveri/ kdt,jdt,inistp,limlow,maxstp,numsum,nummax
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real le, ls          ! latent heat of vaporization/sublimation (j/kg)
      real snow_flux       ! snow heat flux (w/m2)
      real kappa_snow      ! snow thermal conductivity (w/m/k)
      real ice             ! soil layer ice content (mm)
      real ice0            ! last soil layer ice content (mm)
      real delth           ! land layer heat storage change
      real evap            ! evapotranspiration (m/s)
      real ts_old          ! last surface temperature (k)
      real t1_old          ! last first layer soil temperature (k)
      real t1              ! soil temperature (k)
      real tmp_t           ! temporary temperature (k)
      real smmax           ! maximum soil moisture (mm)
      real rad             ! net radiation (w/m2)

      real canopy_evap     ! funcation name (m/s)
      real arno_evap       ! funcation name (m/s)
      real max_unfwat      ! funcation name

      integer m            ! loop index

      real rhoice          ! density of ice (kg/m^3)
      real rhoh2o          ! water density (kg/m^3)
      parameter (rhoice = 917,rhoh2o=1.e3)

      real cp              ! specific heat of moist air at constant (j/kg/k)
      real lf              ! latent heat of freezing (j/kg) at 0c
      real sigma           ! stefan-boltzmann constant (w/m^2/k^-4)
      real t0c             ! ice/water mix temperature (k)
      parameter (cp=1013.0,lf=3.337e5,sigma=5.6730e-8,t0c=273.15)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* initialize latent heat and sensible heat fluxes
* ----------------------------------------------------------------------
      olath = 0.0
      osenh = 0.0
      mfcwt = efcwt
      mfrtf = efrtf

* ----------------------------------------------------------------------
* compute surface temperature at half time step
* ----------------------------------------------------------------------

      if(esnowd .gt. 0.0) then

* ----------------------------------------------------------------------
* compute energy flux through snow pack
* ----------------------------------------------------------------------
         kappa_snow = 2.9302e-6 * efrsn * efrsn
         snow_flux  = kappa_snow * (eftsf-ts) / esnowd
      else
         snow_flux  = 0.0
      endif

* ----------------------------------------------------------------------
* use finite difference method to explicitly solve ground heat
* flux at soil thermal nodes (cherkauer and lettenmaier, 1999)
* ----------------------------------------------------------------------

      ts_old = eftnd(1)
      t1_old = eftnd(2)
      ice0   = efsic(1)

      mftnd(1) = ts             ! the first node of new t_profile is ts

#ifdef DBGVIC
      if(kdt.eq.341)then
      print *,'in func_surfeb check eftnd ', eftnd
      endif
#endif

      call soil_temprfl(nsl,  edtime,   ealpnd,  ebetnd,
     &               egamnd,   eftnd,   ekapnd,   ecsnd,
     &                esmnd,  eicend,   esmxnd,  ebubnd,
     &               eexpnd,   mftnd)

      t1 = mftnd(2)

#ifdef DBGVIC
      if(kdt.eq.341)then
      print *,'in func_surfeb check mftnd ', mftnd
      endif
#endif

* ----------------------------------------------------------------------
* initialize mfsmc
* ----------------------------------------------------------------------

      do m = 1, msl
         mfsmc(m) = efsmc(m)
      end do

* ----------------------------------------------------------------------
* compute the ground heat flux to the top soil layer
* ----------------------------------------------------------------------

      ogrndh = ekapnd(1)*(t1 - ts)/edph(1)
c      if(abs(ogrndh).gt.100) ogrndh = sign(100.0,ogrndh)

* ----------------------------------------------------------------------
* compute the current ice content of the top soil layer
* ----------------------------------------------------------------------

      tmp_t = (ts+ t1)/2.0
      smmax = esmxnd(1)*edph(1)*1000.0

      if(tmp_t .lt. 273.15)then
         ice = efsmc(1) - max_unfwat(tmp_t,smmax,ebubnd(1),eexpnd(1))
         if(ice .lt. 0.) ice=0.
      else
         ice = 0.
      end if

#ifdef DBGVIC
      if(kdt.eq.341)then
      print *,'in func_surfeb ice=',ice, efsmc(1),tmp_t,smmax,
     &         ebubnd(1),eexpnd(1)
      endif
#endif

      delth = ecsnd(1)*((ts_old+t1_old)/2.-tmp_t)*edph(1)/edtime

#ifdef DBGVIC
      if(kdt.eq.341)then
      print *,'delth1=',delth,ecsnd(1),(ts_old+t1_old)/2, tmp_t
      endif
#endif

      delth = delth - rhoice*lf*((ice0-ice)/1000.0)/edtime

* ----------------------------------------------------------------------
* discussion: (ji 2/2004) the contribution from delth to the surface
*    energy balance should be not very high and we should limit the
*    valud of delth. the reason is that in the energy balance equation
*    all the energy terms are facing to a virtual land surface instead
*    of the first soil layer.
* ----------------------------------------------------------------------

      delth = delth/10.0  ! 6/2004
      if(abs(delth).gt.100) delth = sign(100.0, delth)

* temporary set delth=0.0 !ji 2/2004
c      delth = 0.0

#ifdef DBGVIC
      if(kdt.eq.341)then
      print *,'delth2=',delth,' ice=',ice,' ice0=',ice0,efsmc(1)
      endif
#endif

* ----------------------------------------------------------------------
* compute evapotranspiration if not snow covered
* ----------------------------------------------------------------------

      if(esncvfr.lt.1.0) then               ! no snow cover

         rad = enetsht + elongw - sigma*ts**4 + ogrndh + delth

#ifdef DBGVIC
         if(kdt.eq.341)then
         print *,'rad=',rad,enetsht,elongw,sigma*ts**4,ogrndh,delth
         endif
#endif

c         print *,' in func_surfeb efflai',efflai

         if(efflai.gt.0.) then ! snow free canopy

* ----------------------------------------------------------------------
* compute net surface radiation for evaporation estimates
* ----------------------------------------------------------------------

#ifdef DBGVIC
            if(kdt.eq.341)then
            print *,'canopy_evap', msl,  evtype,   edtime,  efflai,
     &                       efcwt,   efrtf,    efsmc,   efsic,
     &                        ewcr,    ewpw,     efrt,     rad,
     &                        evpd,  
     &                       epgcm,   efrtf
            endif
#endif
            evap = canopy_evap(msl,  evtype,   edtime,  efflai,
     &                       efcwt,   efrtf,    efsmc,   efsic,
     &                        ewcr,    ewpw,     efrt,     rad,
     &                     erhoair,  evpair,     evpd,     era,
     &                     enetsht,   etair,    epgcm,   mfcwt,
     &                       mfrtf,   mfsmc)

c            print *,'in func_surfeb canop evap',evap
         else

* ----------------------------------------------------------------------
* compute net surface radiation for evaporation estimates
* ----------------------------------------------------------------------

#ifdef DBGVIC
            if(kdt.eq.341)then
            print *,'arno_evap ',efsmc(1), efsic(1),   smmax, esmr(1),
     &                          ebinf,  edph(1),  edtime,     rad,
     &                         evpair, evpd,
     &		                epgcm
            endif
#endif                     
            mfrtf = efrtf
            evap = arno_evap(      ts,
     &                       efsmc(1), efsic(1),   smmax,  esmr(1),
     &                          ebinf,  edph(1),  edtime,      rad,
     &                        erhoair,   evpair,    evpd,      era,
     &		              enetsht,    etair,   epgcm, mfsmc(1))

c            print *,'in func_surfeb arno evap',evap
         end if

* ----------------------------------------------------------------------
* compute the latent heat flux from the surface and covering vegetation
*   le: latent heat of vaporization (j/kg) 
*   ls: latent heat of sublimation (j/kg)
* ----------------------------------------------------------------------

         if(ts.lt.t0c) then
            ls = (677.-0.07*(ts-t0c))*4.1868*1000.0
            olath = -rhoh2o*ls*evap
         else
            le = (2.501-0.002361*(ts-t0c))*1.e6
            olath = -rhoh2o*le*evap
         endif

#ifdef DBGVIC
         if(kdt.eq.341)then
         print *,'in func_surfeb check',olath,evap
         endif
#endif         

* ----------------------------------------------------------------------
* compute the sensible heat flux from the surface
* ----------------------------------------------------------------------

         osenh = erhoair*cp*(etair - ts)/era
#ifdef DBGVIC
         if(kdt.eq.341)then
         print *,'in func_surfeb osenh ',osenh,erhoair,etair,ts,era,cp
         endif
#endif
      end if

* ----------------------------------------------------------------------
* compute surface energy balance error
* ----------------------------------------------------------------------

      func_surfeb = (1.0-esncvfr)*(enetsht + elongw - sigma*ts**4 +
     &              olath + osenh) + ogrndh + delth +
     &              esncvfr*snow_flux


#ifdef DBGVIC
      if(kdt.gt.339)then
      print *,'in func_surfeb esnowd(m) =',esnowd,esncvfr,enetsht,
     &       elongw,sigma*ts**4,ts,olath,osenh,ogrndh,delth,snow_flux
      print *,'in func_surfeb func_surfeb=',func_surfeb
      endif
#endif

#endif

      return
      end
