      function max_unfwat(t, soilmx, bubble, expt)

#include <vartyp.h>

************************************************************************
* source file:       max_unfwat.f
* purpose:           compute the maximum amount of unfrozen 
*                    water that can exist at the current temperature.
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic c version
*                    (maximum_unfrozen_water in soil_conduction.c)
* notes:
*     please check equnation (14) in cherkauer and lettenmaier 
*     (jgr 1999) for the equation of computing unfrozen soil moisture
************************************************************************

* ------------------- input variables ----------------------------------
      real t          ! soil temperature (k)
      real soilmx     ! maximum soil moisture (mm or m3/m3)
      real bubble     ! bubbling pressure of soil (cm)
      real expt       ! parameter for ksat with soil moisture (n/a)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real max_unfwat ! unfrozen soil moisture (mm or m3/m3)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real unfrozen, bas

      real g          ! gravity acceleration (m/s2)
      real lf         ! latent heat of freezing (j/kg) at 0c
      parameter (g=9.81, lf = 3.337e5)
* ----------------------------------------------------------------------

      bas = (-lf * (t-273.15)) / t  / (g * bubble / 100.)
      unfrozen = soilmx * bas **( -2.0/(expt-3.0))

      if(unfrozen .gt. soilmx) unfrozen = soilmx
      if(unfrozen .lt. 0) unfrozen = 0
      
      max_unfwat = unfrozen

      return
      end
