      function stabilityc(z, d, tsurf, tair, wind, z0)

#include <vartyp.h>

************************************************************************
* source file:       stabilityc.f
* purpose:           this function calculates the stability correction for 
*                    exchange of sensible heat between surface&atmosphere
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (stabilitycorrection.c)
************************************************************************

* ------------------- input variables ----------------------------------
      real z          ! reference height (m)
      real d          ! displacement height (m)
      real tsurf      ! surface temperature (k)
      real tair       ! air temperature (k)
      real wind       ! wind speed (m/s)
      real z0         ! roughness length (m)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real stabilityc ! multiplier for aerodynamic resistance
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real rilimit    ! upper limit for richardsons number
      real ri         ! richardsons number 
      real ricr       ! critical richardsons number
      real g          ! gravity acceleration (m/s2)

      parameter (ricr = 0.2, g=9.81)
*-----------------------------------------------------------------------

* ----------------------------------------------------------------------
* calculate effect of atmospheric stability using richardson number approach
* ----------------------------------------------------------------------

      stabilityc = 1.0

      if (tsurf.ne.tair.and.z.gt.d) then

* ----------------------------------------------------------------------
* non-neutral conditions
* ----------------------------------------------------------------------

         ri = g*(tair-tsurf)*(z-d)/((tair+tsurf)/2.0*wind*wind)
    
         if(z0.le.0.0) z0=0.01

         if(((z-d)/z0).le.1) then
            rilimit = tair/((tair+tsurf)/2.0 * 5.0)
         else
            rilimit = tair/((tair+tsurf)/2.0*(log((z-d)/z0) + 5.0))
         end if

         if (ri .gt. rilimit) ri = rilimit
    
         if (ri .gt. 0.0) then 
            stabilityc = (1 - ri/ricr) * (1 - ri/ricr)
         else      
            if (ri .lt. -0.5) ri = -0.5
            stabilityc = sqrt(1.0 - 16.0 * ri)
         end if
      end if

#ifdef DBGVIC
      print *,'vic in stability check stability ', stabilityc
#endif

      return
      end
