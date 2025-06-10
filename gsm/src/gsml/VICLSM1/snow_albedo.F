      function snow_albedo(snowage, tsf)

#include <vartyp.h>

************************************************************************
* source file:       snow_albedo.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (snow_utility.c)
*
* abstract: this function computes the snow pack surface albedo based 
*           on snow age and season, using the tables generated in 
*           snow_table_albedo.
************************************************************************

* ------------------- input variables ----------------------------------
      real snowage     ! surface snow age (day)
      real tsf         ! temperature of snow pack surface layer (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real snow_albedo ! snow albedo
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real t0c
      parameter (t0c = 273.15)

      real new_snow_alb
      parameter (new_snow_alb = 0.85)

      real snow_alb_accum_a, snow_alb_accum_b
      parameter (snow_alb_accum_a=0.94, snow_alb_accum_b=0.58)

      real snow_alb_thaw_a, snow_alb_thaw_b
      parameter (snow_alb_thaw_a=0.82, snow_alb_thaw_b=0.46)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* new snow 
* ----------------------------------------------------------------------
      if(snowage .eq. 0.0) then
         snow_albedo = new_snow_alb

* ----------------------------------------------------------------------
* aged snow: accumulation season 
* ----------------------------------------------------------------------
      else
         if(tsf .lt. t0c) then
            snow_albedo = new_snow_alb*snow_alb_accum_a**
     &                    (snowage**snow_alb_accum_b)

* ----------------------------------------------------------------------
* melt season
* ----------------------------------------------------------------------
         else
            snow_albedo = new_snow_alb*snow_alb_thaw_a** 
     &                    (snowage**snow_alb_thaw_b)
         end if
      end if

      return
      end
