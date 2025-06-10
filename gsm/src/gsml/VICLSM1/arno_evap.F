      function arno_evap(   ts,
     &                   fsmc1,  soili,  smmax,  smres,
     &                    binf,     d1,  dtime,    rad,
     &                  rhoair,  vpair,    vpd,     ra,
     &                  netsht,   tair,   pgcm, mfsmc1)

#include <vartyp.h>

************************************************************************
* source file:       arno_evap.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (arno_evap.c)
*
* abstract: this routine compute evaporation based on the assumption that 
*    evaporation is at the potential for the area which is saturated, and at
*    some percentage of the potential for the area which is partial saturated.
************************************************************************

* ------------------- input variables ----------------------------------
      real ts             ! surface temperature (k)
      real soili          ! soil ice in layer 1 (mm)
      real smmax          ! maximum soil moisture (mm)
      real smres          ! residual soil moisture (mm)
      real binf           ! vic infiltration parameter (n/a)
      real d1             ! soil layer 1 thickness (m)
      real dtime          ! time step (second)
      real rad            ! net radiation (w/m2)
      real rhoair         ! air density (kg/m3)
      real vpair          ! actual vapor pressure of air (pa)
      real vpd            ! vapor pressure deficit (pa)
      real ra             ! aerodynamical resistance (s/m)
      real netsht         ! net shortwave radiation (w/m2)
      real tair           ! air temperature (k)
      real pgcm           ! pressure (pa)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fsmc1          ! soil moisture in layer 1 (mm)
      real mfsmc1         ! modified soil moisture in layer 1 (mm)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real arno_evap      ! evaporation over bare soil (m/s)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real evap           ! evaporation
      real moist          ! available soil moisture (mm)
      real epot           ! potential bare soil evaporation (m/s)
      real max_infil      ! maximum infitration capacity (mm)
      real tmp            ! infiltration capacity with moist (mm)
      real tmp_evap       ! temporary varible for evaporation (mm/time step)
      real ratio          ! for computing tmp
      real as             ! saturated area (fraction)
      real beta_asp
      real dummy          ! temporary variable
      real tmpsum

      integer num_term, i ! loop index

      real vicpenman, svp ! function name
      real rhoh2o         ! water density (kg/m^3)
      parameter (rhoh2o=1.e3)
* ----------------------------------------------------------------------

      evap = 0

      moist = fsmc1 - soili
      if(moist.gt.smmax) moist= smmax

* ----------------------------------------------------------------------
* calculate the potential bare soil evaporation (m/s)
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* in the original c-vic (uw), the negative rad is also used in penman
* method, which will introduce errors in computing energy balance.
* therefore, we control the negatvie rad computation. if rad is negative,
* we use aerodynamic method to compute evap and assume canopy 
* temperature equals to air temperature.
* ----------------------------------------------------------------------

      if(rad.gt.0.0)then
         epot = vicpenman(rad,     vpd,      ra,      0.0, 
     &                    0.0,     1.0,     1.0,     tair, 
     &                 netsht,    pgcm,       0)
      else
         epot = rhoair*(0.622/pgcm)*(svp(ts)-vpair)/ra
         epot = epot / rhoh2o
         if (vpd.eq.0.0.and.epot.gt.0.0) epot=0.0
      endif

#ifdef DBGVIC
      print *,'in arno_evap epot=',epot
#endif

* ----------------------------------------------------------------------
* compute temporary infiltration rate based on given soil_moist.
* ----------------------------------------------------------------------

      max_infil = (1.0+binf)*smmax
      if(binf .eq. (-1.0) .or. moist.eq.smmax) then
         tmp = max_infil
      else 
         ratio = 1.0 - moist /smmax

* ----------------------------------------------------------------------
* if(ratio < small && ratio > -small) ratio = 0.
* ----------------------------------------------------------------------
         ratio = ratio**(1.0 / (binf + 1.0))
         tmp = max_infil*(1.0 - ratio)       !infiltration capacity
      endif

* ----------------------------------------------------------------------
* evaporation see eqs.(14) & (15) in liang et al. (jgr 1994) derivation.
* ----------------------------------------------------------------------

      if(tmp .ge. max_infil.or.epot.le.0.0) then
         evap = epot
      else

* ----------------------------------------------------------------------
* compute as. as is % area saturated, 1-as is % area that is unsaturated
* ----------------------------------------------------------------------
         ratio = tmp/max_infil 
         ratio = 1.0 - ratio
         ratio = ratio**binf
         as = 1 - ratio                 ! saturated area

* ----------------------------------------------------------------------
* compute the beta function in the arno evaporation model using
* the first 30 terms in the power expansion expression.
* ----------------------------------------------------------------------

         ratio = ratio**(1.0/binf)
         dummy = 1.0
         do num_term = 1, 30
            tmpsum = ratio
            do i = 1, num_term
               tmpsum = tmpsum*ratio
            end do
            dummy = dummy + binf*tmpsum/(binf+num_term)
         end do
            
         beta_asp = as + (1.0-as)*(1.0-ratio)*dummy
         evap = epot*beta_asp
      endif

#ifdef DBGVIC
      print *,'in arno_evap evap=',evap
#endif

* ----------------------------------------------------------------------
* evaporation cannot exceed available soil moisture.
* evaporation second soil layer = 0.0
* ----------------------------------------------------------------------

      tmp_evap = evap*1000.0*dtime   ! m/s -> mm/time step

      if((moist-smres).gt.0.0) then
         if(tmp_evap .gt. (moist - smres)) then
            evap = (moist -  smres)/1000.0/dtime
            tmp_evap=moist-smres
         end if
      else
         tmp_evap = 0.0
         evap = 0.0
      endif

#ifdef DBGVIC
      print *,'in arno_evap tmp_evap',tmp_evap,moist,smres,evap
#endif

      mfsmc1 = fsmc1 - tmp_evap
c      mfsmc1 = fsmc1
      arno_evap = evap       ! m/s

      return
      end
