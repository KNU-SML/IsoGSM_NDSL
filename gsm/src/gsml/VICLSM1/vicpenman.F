      function vicpenman(rad,     vpd,      ra,      rs, 
     &                  rarc,     lai, gsm_inv,    tair, 
     &              netshort,    pgcm,     rgl)

#include <vartyp.h>

************************************************************************
* source file:       vicpenman.f
* purpose:           calculate et using the combination eq
* prgmmr:            ji chen
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic c version (penman.c)
*
* comments:
*    please check equations (11), (12), (13), (14), (15), (16), and (3)
*    in wigmosta, vail and lettenmaier (wrr 1994) for computing evap.
************************************************************************

* ------------------- input variables ----------------------------------
      real rad        ! net radiation (w/m2)
      real vpd        ! vapor pressure deficit (pa)
      real ra         ! aerodynamical resistance (s/m)
      real rs         ! minimum stomatal resistance (s/m)
      real rarc       ! architectural resistance (s/m)
      real lai        ! lai
      real gsm_inv    ! soil moisture stress factor
      real netshort   ! net shortwave radiation (w/m2)
      real tair       ! air temperature (k)
      real pgcm       ! pressure (pa)
      real rgl        ! value of solar radiation below which there 
                      ! will be no transpiration (ranges from ~30 w/m^2
                      ! for trees to ~100 w/m^2 for crops)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real vicpenman  ! evapotranspiration (m/s)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real ta         ! air temperature (c)
      real f
      real tmp_value
      real slope      ! slope of saturated vapor pressure curve (pa/c)
      real rc         ! canopy resistance
      real rhoair     ! density of air in (kg/m3)
      real le         ! latent heat of vaporization (j/kg)
      real gamma      ! psychrometric constant (pa/c)
      real tfactor    ! factor for canopy resist based on temperature
      real vpdfactor  ! factor for canopy resistance based on vpd
      real dayfactor  ! factor for canopy resist based on photosynthesis
      real svp_slope  ! function name

      real cp         ! specific heat of moist air at constant (j/kg/k)
      real t0c        ! ice/water mix temperature (k)
      parameter (cp=1013.0, t0c=273.15)

      real closure    ! pa
      real rsmax      ! maximum stomatal resistance (s/m)
      real vpdminftr  ! minimum vpd factor
      parameter (closure = 4000.0,rsmax=5000,vpdminftr=0.1)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* calculate the slope of the saturated vapor pressure curve in pa/k
* ----------------------------------------------------------------------

      slope = svp_slope(tair)

      ta = tair - t0c

* ----------------------------------------------------------------------
* calculate resistance factors (wigmosta et al., 1994)
* ----------------------------------------------------------------------

      if(rs .gt. 0.) then
         f = netshort / rgl
         dayfactor = (1. + f)/(f + rs/rsmax)
      else
         dayfactor = 1.
      end if

      tfactor = .08 * ta - 0.0016 * ta * ta
      if(tfactor.le.0.0) tfactor = 1.e-10

      vpdfactor = 1 - vpd/closure
      if(vpdfactor .lt. vpdminftr) vpdfactor = vpdminftr

* ----------------------------------------------------------------------
* calculate canopy resistance in s/m
* ----------------------------------------------------------------------

      tmp_value = lai * gsm_inv * tfactor * vpdfactor

      if(tmp_value.ne.0.0) then
         rc = rs/(lai * gsm_inv * tfactor * vpdfactor) * dayfactor
         if(rc .gt. rsmax) rc = rsmax
      else
         rc = rsmax
      end if

c      print *,'penman rc',rc,rs,lai,gsm_inv
c      print *,'penman tfactor ',tfactor,ta,vpdfactor,vpd
c      print *,'penman dayfactor',dayfactor,netshort,rgl

* ----------------------------------------------------------------------
* calculate latent heat of vaporization. eq. 4.2.1 in handbook of 
* hydrology, assume ts is ta.  le (j/kg)
* ----------------------------------------------------------------------
      le = (2.501 - 0.002361 * ta) * 1.0e6 
  
* ----------------------------------------------------------------------
* calculate gamma (pa/c) eq. 4.2.28. handbook of hydrology 
* ----------------------------------------------------------------------
      gamma = 1628.6 * pgcm/le
  
* ----------------------------------------------------------------------
* calculate air density (kg/m3), using eq. 4.2.4 handbook of hydrology 
* ----------------------------------------------------------------------
      rhoair = 0.003486 * pgcm/(275.0 + ta)
 
* ----------------------------------------------------------------------
* calculate the evaporation in mm/s (not dividing by the density 
* of water (~1000 kg/m3)), the result ends up being in mm instead of m 
* ----------------------------------------------------------------------
    
      vicpenman = (slope * rad + rhoair * cp * vpd/ra)/
     &    (le * (slope + gamma * (1 + 0.1*(rc + rarc)/ra)))

c      print *,'vpenman=',vicpenman
c      print *,'vicpenman ',slope,rad,rhoair,cp,vpd,ra,rc,rarc,gamma,le
c      print *,'penman ',slope * rad,rhoair * cp * vpd/ra,
c     &                  gamma * (1 + (rc + rarc)/ra)

      vicpenman = vicpenman / 1000.0         ! mm/s -> m/s
      if (vpd .gt. 0.0 .and. vicpenman .lt. 0.0) vicpenman = 0.0 ! no dew

      return
      end 
