      function heat_capcty(soil_fract, water_fract, ice_fract)

#include <vartyp.h>

************************************************************************
* source file:       heat_capcty.f
* purpose:           this function calculates the soil volumetric heat 
*                    capacity based on the fractional volume of its 
*                    component parts.
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              november 2003
* prgm history:      modified from 4.0.3 uw vic c version
*                    (volumetric_heat_capacity in soil_conduction.c)
* notes:
*       constant values are volumetric heat capacities in j/m^3/k
*       soil value is for clay or quartz - assumed for all other types
************************************************************************

* ------------------- input variables ----------------------------------

      real soil_fract    ! fraction of soilvolume composed of actual soil
      real water_fract   ! fraction of soilvolume composed of liquid water
      real ice_fract     ! fraction of soilvolume composed of ice
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real heat_capcty   ! heat capacity (j m-3 k-1)
*-----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real cs
      real organic_fract
*-----------------------------------------------------------------------

      organic_fract = 0.0

      cs = 2.0e6 * (soil_fract - organic_fract)
      cs = cs + 4.2e6 * water_fract
      cs = cs + 1.9e6 * ice_fract
      cs = cs + 2.7e6 * organic_fract
      cs = cs + 1.3e3 * (1.-(soil_fract+water_fract+
     &     ice_fract+organic_fract))

      heat_capcty = cs

      return
      end 
