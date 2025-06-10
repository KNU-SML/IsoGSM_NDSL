      function svp(temp)

#include <vartyp.h>

************************************************************************
* source file:       svp.f
* purpose:           compute the saturated vapor pressure using 
*                    handbook of hydrology eqn 4.2.2 pressure in pa
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (svp.c)
************************************************************************

* --------------- input variables from atmospheric model ---------------
      real temp        ! air temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real svp         ! saturated vapor pressure (pa)
* ----------------------------------------------------------------------

* --------------- local variables --------------------------------------
      real t           ! air temperature (c)
      real a_svp, b_svp, c_svp, t0c
      parameter (a_svp=0.61078, b_svp=17.269, c_svp=237.3, t0c=273.15)
* ----------------------------------------------------------------------

      t = temp - t0c

      svp = a_svp * exp((b_svp * t)/(c_svp+t))

      if(t.lt.0) then
         svp = svp*(1.0 + .00972 * t + .000042 * t * t)
      end if

      svp = svp * 1000.0

      return
      end
