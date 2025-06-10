      function snow_density(new_snow, tair, swq, depth, dtime, tsf)

#include <vartyp.h>

************************************************************************
* source file:       snow_density.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august, & september 2003
* prgm history:      modified from 4.0.3 uw vic (in snow_utility.c)
*
* abstract: this function computes the snow density based on the day of 
*           the year. density information comes from a plot of seasonal 
*           variation of typical snow densities found in bras 
*           (figure 6.10, p 258). the equation was developed by regressing 
*           against the curve for southern manitoba, so this routine should 
*           be modified if used outside the plains of south central canada, 
*           and the north central us. 
************************************************************************

* ------------------- input variables ----------------------------------
      real dtime      ! model time step
      real new_snow   ! amount of new snow (m/time step)
      real tair       ! air temperature (k)
      real swq        ! snow water equivalent (m)
      real depth      ! snow depth (m)
      real tsf        ! snow surface temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real snow_density ! (kg/m3)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real new_sdens       ! new snow density(kg/m3)
      real max_change      ! maximum fraction of snowpack depth change
      parameter (new_sdens=50., max_change=0.9)

      real eta0            ! viscosity of snow at t = 0c and density = 0
                           ! used in calculation of true viscosity (ns/m2)
      real c5              ! constant used in snow viscosity calculation,
			   ! taken from snthrm.89 (/c)
      real c6              ! constant used in snow viscosity calculation
      real g               ! gravitational accelleration (m/(s^2))
      parameter (eta0 = 3.6e6, c5 = 0.08, c6 = 0.021, g=9.81)

      real t0c, rhoh2o
      parameter (t0c=273.15, rhoh2o=1.e3)

      real ta, depth_new
      real density, pdepth, density_new
      real overburden, viscosity, ddepth
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* compaction of snow pack by new snow fall *** bras pg. 257 *
* ----------------------------------------------------------------------

      if(new_snow .gt. 0) then

* ----------------------------------------------------------------------
* estimate density of new snow based on air temperature 
* ----------------------------------------------------------------------

         ta = (tair-t0c) * 9. / 5. + 32.
         if(ta .gt. 0) then
            density_new = new_sdens + 1000.
     &               * (ta / 100.) * (ta / 100.)
         else
            density_new = new_sdens
         end if

         if(depth .gt. 0.) then

* ----------------------------------------------------------------------
* compact current snowpack by weight of new snowfall
* ----------------------------------------------------------------------

           pdepth = ((new_snow / 0.0254)*(depth / 0.0254)) / 
     &          (swq/0.0254)* (((depth/0.0254)/10.)**0.35) * 0.0254
  
* ----------------------------------------------------------------------
* pdepth cannot be greater than depth 
* ----------------------------------------------------------------------

            if (pdepth .ge. depth) then
               pdepth = max_change * depth
            end if

            depth_new = 1000.0 * new_snow / density_new
            depth = depth - pdepth + depth_new
            swq   = swq + new_snow
            if(depth.gt.0) then
               density = 1000. * swq / depth
            else
               density = density_new
            endif
         else

* ----------------------------------------------------------------------
* no snowpack present, so snow density equals that of new snow
* ----------------------------------------------------------------------
            density = density_new
            swq = swq + new_snow
         end if
      else
         if(depth.gt.0) then
            density = 1000. * swq / depth
         else
            density = density_new
         endif
      end if

* ----------------------------------------------------------------------
* densification of the snow pack due to aging
* based on snthrm89 r. jordan 1991 - used in barts dhsvm code
* ----------------------------------------------------------------------

      if(swq.gt.0) then
         depth       = 1000. * swq / density
         overburden  = 0.5 * g * rhoh2o * swq
         viscosity   = eta0 * exp(-c5 * (tsf-t0c) + c6 * density)
         ddepth      = -overburden / viscosity * depth * dtime
         if(abs(ddepth).gt.0.05*depth*dtime/86400.0) then
            ddepth = - 0.05*depth*dtime/86400.0
         endif
         depth       = depth + ddepth
         
         if(depth .le. 0) depth = (1-max_change)*depth
         
         density     = 1000. * swq / depth
      endif

      snow_density = density

      return
      end
