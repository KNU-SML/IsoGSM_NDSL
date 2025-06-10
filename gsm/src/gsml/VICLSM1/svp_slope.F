      function svp_slope(temp)

#include <vartyp.h>

************************************************************************
* source file:       svp_slope.f
* purpose:           the gradient of d(svp)/dt using handbook
*                    of hydrology eqn 4.2.3
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (svp.c)
************************************************************************

* --------------- input variables from atmospheric model ---------------
      real temp        ! air temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real svp_slope   ! the gradient of saturated vapor pressure (pa/c)
* ----------------------------------------------------------------------

* --------------- local variables --------------------------------------
      real t           ! air temperature (c)
      real svp         ! function name
      real b_svp, c_svp, t0c
      parameter (b_svp = 17.269, c_svp = 237.3, t0c=273.15)
* ----------------------------------------------------------------------

      t = temp - t0c
      
      svp_slope = (b_svp * c_svp) / ((c_svp + t) * (c_svp + t)) 
     &            * svp(temp)

      return
      end
