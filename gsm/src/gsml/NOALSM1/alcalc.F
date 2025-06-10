      subroutine alcalc (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
c
#include <define.h>
      implicit none
c      
c ----------------------------------------------------------------------
c calculate albedo including snow effect (0 -> 1)
c   alb     snowfree albedo
c   snoalb  maximum (deep) snow albedo
c   shdfac    areal fractional coverage of green vegetation
c   shdmin    minimum areal fractional coverage of green vegetation
c   sncovr  fractional snow cover
c   albedo  surface albedo including snow effect
c   tsnow   snow surface temperature (k)
c ----------------------------------------------------------------------
      real alb, snoalb, shdfac, shdmin, sncovr, albedo, tsnow
c ----------------------------------------------------------------------
c snoalb is argument representing maximum albedo over deep snow,
c as passed into sflx, and adapted from the satellite-based maximum 
c snow albedo fields provided by d. robinson and g. kukla 
c (1985, jcam, vol 24, 402-411)
c ----------------------------------------------------------------------
#ifndef NOAH_JUN05
          albedo = alb + (1.0-(shdfac-shdmin))*sncovr*(snoalb-alb) 
#else
          albedo=alb+sncovr*(snoalb-alb)
#endif
          if (albedo .gt. snoalb) albedo=snoalb
c
c     base formulation (dickinson et al., 1986, cogley et al., 1990)
c          if (tsnow.le.263.16) then
c            albedo=snoalb
c          else
c            if (tsnow.lt.273.16) then
c              tm=0.1*(tsnow-263.16)
c              albedo=0.5*((0.9-0.2*(tm**3))+(0.8-0.16*(tm**3)))
c            else
c              albedo=0.67
c            endif
c          endif
c
c     isba formulation (verseghy, 1991; baker et al., 1990)
c          if (tsnow.lt.273.16) then
c            albedo=snoalb-0.008*dt/86400
c          else
c            albedo=(snoalb-0.5)*exp(-0.24*dt/86400)+0.5
c          endif
c
c ----------------------------------------------------------------------
c end subroutine alcalc
c ----------------------------------------------------------------------
      return
      end
