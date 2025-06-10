      subroutine massrelease(csnow, tmpintstr, releasmss, drip)

#include <vartyp.h>

************************************************************************
* source file:       massrelease.f
* purpose:           calculates mass release of snow from canopy
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic 
*                    (massrelease.c and part of dhsvm)
************************************************************************

* ------------------- input/output variables ---------------------------
      real csnow         ! canopy intercepted snow (m)
      real tmpintstr     ! temporary intercepted storage (m)
      real releasmss     ! released snow mass from leaf (m)
      real drip          ! leaf drip (m)

* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real min_is        ! min_interception_storage (m)
      parameter (min_is = 0.005)

      real threshold, maxrelease, tempdrip, tmprelmass
*-----------------------------------------------------------------------
  
* ----------------------------------------------------------------------
* if the amount of snow in the canopy is greater than some minimum
* value, min_is, then calculate mass release and drip
* ----------------------------------------------------------------------

 5000 if (csnow .gt. min_is) then
         threshold  = 0.10 * csnow
         maxrelease = 0.17 * csnow
    
* ----------------------------------------------------------------------
* if the amount of snow_melt after interception, snow_melt, is >= the
*  theshhold then there is mass release.  if snow_melt is < the treshhold
*  then there is no mass release but that water remains in
*  tempintstorage which will be augmented during the next
*  compute period
* ----------------------------------------------------------------------
    
         if (tmpintstr .ge. threshold) then
            drip  = drip  + threshold
            csnow = csnow - threshold
            tmpintstr = tmpintstr - threshold

            if (csnow .lt. min_is) then
               tmprelmass = 0.0
            else
               tmprelmass = min((csnow-min_is),maxrelease)
            end if

            releasmss = releasmss + tmprelmass
            csnow = csnow - tmprelmass
            go to 5000
         else
            tempdrip = min(tmpintstr, csnow)
            drip  = drip + tempdrip
            csnow = csnow - tempdrip
         end if
      else

* ----------------------------------------------------------------------
* interceptedsnow < min_is) if the amount of snow in
*   the canopy is less than some minimum value, min_is,
*   then only melt can occur and there is no mass release.
* ----------------------------------------------------------------------

         tempdrip = min(tmpintstr, csnow)
         drip = drip + tempdrip
         csnow = csnow - tempdrip
         tmpintstr = 0.0
      end if

      return
      end
