      subroutine snow_intcept(dtime,     lai, rhoair,    tair,
     &                        vpair,     vpd,     ra,    wind,
     &                       netsht,   longw,  rainf,   snowf,
     &                         pgcm,  cwater,  csnow,  raintf,
     &                       snowtf,  lathsc, senhsc)

#include <vartyp.h>

************************************************************************
* source file:       snow_intercept.f
* purpose:           calculates the interception and subsequent release
*                    of snow by the forest canopy using energy balance
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (snow_intercept.c)
*
* comments: only the top canopy layer is taken into account for snow
*           interception. snow interception by lower canopy is
*           disregarded. rain water can be intercepted by lower canopy
*           layers. of course:  no vegetation -> no interception
************************************************************************

* ------------------- input variables ----------------------------------
      real dtime         ! time step (second)
      real lai           ! vegetation leaf area index
      real rhoair        ! air density (kg/m^3)
      real tair          ! air temperature (k)
      real vpair         ! actual vapor pressure of air (pa)
      real vpd           ! air vapor pressure deficit (pa)
      real ra            ! aerodynamic resistance (s/m)
      real wind          ! wind speed related to ra(2) (m/s)
      real netsht        ! net short wave (w/m^2)
      real longw         ! long wave (w/m^2)
      real rainf         ! rain fall (m/time step)
      real snowf         ! snow fall (m/time step)
      real pgcm          ! air pressure (pa)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real cwater        ! canopy intercepted water (m)
      real csnow         ! canopy intercepted snow (m)
* ----------------------------------------------------------------------

* ---------------------- output variables ------------------------------
      real raintf        ! rain throughfall (m/time step)
      real snowtf        ! snow throughfall (m/time step)
      real lathsc        ! leaf snow latent heat (w/m2) (+: air to snow)
      real senhsc        ! leaf snow sensible heat flux (w/m2) (+: to snow)
*-----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real imax1         ! maxium water intecept regardless of temp (m)
      real maxint        ! maximum rainf interception storage (m)
      real blownsnow     ! depth of snow blown of the canopy (m) 
      real delsnwint     ! change in swe of snow interceped on branches (m)
      real maxwatint     ! water interception capacity (m)  
      real maxsnwint     ! snow interception capacity (m) 
      real intrf         ! fraction of intercepted water which is liquid
      real intsf         ! fraction of intercepted water which is solid
      real overload      ! temp variable to calculated structural overloading

      real ta            ! air temperature (c)
      real tcanopy       ! canopy temperature (k)
      real longout       ! emitted longwave radiation (w/m2)
      real advectede     ! advaected energy (w/m2)
      real rad           ! net radiation (w/m2)
      real essnow        ! vapor pressure (pa)
      real vaporflux     ! vapor flux from canopy (m/s) [+ air -> canpoy]
      real vaporf        ! vapor from canopy (m) (+ air -> canopy)
      real ls            ! latent heat of sublimation (j/kg)
      real refreezee     ! refreeze energy (j)
      real potsnwmt      ! potential snow melt (m)
      real exssnwmt      ! excess snow melt (m)
      real tmpintstr     ! temporary intercepted storage (m)
      real releasmss     ! released snow mass from leaf (m)
      real drip          ! leaf drip (m)
      real inicsnow      ! initial canopy snow storage (m)
      real tmp_csw       ! temporary canopy water (m)
      real svp           ! function name

      real cpair         ! specific heat at constant pressure of air (j/deg/k)
      parameter (cpair=1004.0)

      real lai_sm        ! lai snow multiplier (m)
      real lai_wf        ! leaf water factor for interception storage (m)
      real liq_wc        ! liquid water capacity
      real min_is        ! min snow interception storage (m)
      parameter (lai_sm=0.0005,lai_wf=0.002,liq_wc=0.035,min_is=0.005)

      real t0c           ! ice/water mix temperature (k)
      real sigma         ! stefan-boltzmann constant (w/m^2/k^-4)
      real ch_water      ! volumetric heat capacity (j/(m3*c)) of water
      parameter (t0c=273.15, sigma=5.67e-8,ch_water = 4186.8e3)

      real lf            ! latent heat of freezing (j/kg) at 0c
      real rhoh2o        ! water density (kg/m^3)
      parameter (lf=3.337e5,rhoh2o=1.e3)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* initialize variables
* ----------------------------------------------------------------------
      raintf = 0.0
      snowtf = 0.0

      imax1  = 4.0* lai_sm * lai          ! m
      maxint = lai_wf * lai               ! m

* ----------------------------------------------------------------------
* during snow falling, all the snow and water intercepted by vegetation,
* and soil moisture are averaged over the whole grid cell,
* and we only compute the interception over the wet district (mu=1.0);
* but we still need separate canopy snow interception into the wet
* and dry districts in case of rain fall after snowing, resulting in 
* heterogeneity of snow coverage. ji 2003
* ----------------------------------------------------------------------
         
      if(snowf .gt. 0.) then  ! calculate snow interception

* ----------------------------------------------------------------------
* determine the maximum snow interception water equivalent.           
* kobayashi, d., 1987, snow accumulation on a narrow board,           
* cold regions science and technology, (13), pp. 239-245. figure 4.
* ----------------------------------------------------------------------

         ta = tair - t0c

         if (ta.lt.-1.0.and.ta.gt.-3.0)then
            maxsnwint = (ta*3.0/2.0) + (11.0/2.0)
         else if (ta.gt.-1.0) then
            maxsnwint = 4.0
         else
            maxsnwint = 1.0
         end if
         
* ----------------------------------------------------------------------
* therefore lai_ratio decreases as temp decreases
* ----------------------------------------------------------------------

         maxsnwint = maxsnwint*lai_sm*lai
  
* ----------------------------------------------------------------------
* calculate snow interception
* ----------------------------------------------------------------------

         delsnwint = (1.0-csnow/maxsnwint)*snowf
         if ((delsnwint+csnow).gt.maxsnwint)then
            delsnwint = maxsnwint - csnow
         end if
         if (delsnwint .lt. 0.0) delsnwint = 0.0
            
* ----------------------------------------------------------------------
* reduce the amount of intercepted snow if snowing, windy and cold.
* (< -3 to -5 c). schmidt and troendle 1992 western snow conference paper.
* ringyo shikenjo tokyo, #54, 1952. bulletin of the govt. 
* forest exp. station, govt. forest exp. station, meguro, tokyo, japan.
* forstx 634.9072 r475r #54. page 146, figure 10. 
* ----------------------------------------------------------------------

         if(ta.lt.-3.0.and.delsnwint.gt.0.0.and. wind.gt.1.0) then
            blownsnow = (0.2 * wind-0.2)*delsnwint
            if(blownsnow.ge.delsnwint) blownsnow = delsnwint
            delsnwint = delsnwint - blownsnow
         end if
  
* ----------------------------------------------------------------------
* now update snowf and total accumulated intercepted snow amounts
* ----------------------------------------------------------------------

         if ((csnow+delsnwint).gt.imax1) delsnwint=0.0 
  
* ----------------------------------------------------------------------
* pixel depth snow through fall (snowtf)
* ----------------------------------------------------------------------

         snowtf = snowf - delsnwint

* ----------------------------------------------------------------------
* physical depth
* ----------------------------------------------------------------------
         csnow  = csnow + delsnwint
      end if

      if(rainf.gt.0) then

* ----------------------------------------------------------------------
* calc amount of rain intercepted on branches and stored in intercepted snow
*
* if snow fall and rain fall mixing, newmu = 1.0, just compute rain 
* interception over the wet district; otherwise, we need both the 
* dry and wet districts. however, since there is now rain fall over
* the dry district, we also only compute rain fall interception over
* the wet district. ji 2003
*
* before a new rain fall event occuring, the canopy water and rainthrough 
* fall are averaged over the whole grid cell; during the rain fall event, 
* the canopy water and rainthrough fall are updated according to the 
* fraction of the wet district on where rain fall falls. ji 2003
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* physical depth
* ----------------------------------------------------------------------

         maxwatint=liq_wc*csnow+maxint
            
         if ((cwater+rainf).le.maxwatint)then
            cwater = cwater + rainf                ! physical depth
         else
            raintf = cwater + rainf - maxwatint
            cwater = maxwatint                     ! physical depth
         end if
      end if

* ----------------------------------------------------------------------
* at this point we have calculated the amount of snow fall intercepted
* and the amount of rain fall intercepted. these values have been 
* appropriately subtracted from snow fall and rain fall to determine 
* snow throughfall and rain throughfall. however, we can end up with the 
* condition that the total intercepted rain plus intercepted snow is 
* greater than the maximum bearing capacity of the tree regardless of air 
* temp (imax1). the following routine will adjust intrain and intsnow 
* by triggering mass release due to overloading. of course since intrain
* and intsnow are mixed, we need to slough them off as fixed fractions  
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* trigger structural unloading
* ----------------------------------------------------------------------

      if(lai.gt.0.0) then
         tmp_csw = cwater + csnow
         if (tmp_csw.gt.imax1)then 
            overload = tmp_csw - imax1
            intrf = cwater/tmp_csw
            intsf = csnow/tmp_csw
            cwater= cwater - overload*intrf
            csnow = csnow  - overload*intsf
            raintf= raintf + overload*intrf
            snowtf= snowtf + overload*intsf
         end if
      end if

      inicsnow = csnow

* ----------------------------------------------------------------------
* the canopy temperature is assumed to be equal to the air temperature if 
* the air temperature is below t0c, otherwise the canopy temperature is 
* equal to t0c (273.15 k)
* ----------------------------------------------------------------------

      if(tair .gt. t0c) then
         tcanopy = t0c
      else
         tcanopy = tair
      end if

* ----------------------------------------------------------------------
* calculate the net radiation at the canopy surface, using the canopy 
* temperature.  the outgoing longw is subtracted twice, because the 
* canopy radiates in two directions
* ----------------------------------------------------------------------

      longout = sigma * tcanopy**4
      rad     = netsht + longw - 2*longout
      if(rad.lt.0.) then        ! in case of setting too high tcanopy
         rad = 0.0
         tcanopy = ((netsht + longw)/2.0/sigma)**(0.25)
         if(tcanopy.lt.(min(t0c,tair)-5.0)) then
            tcanopy = min(t0c,tair)-5.0
            longout = sigma * tcanopy**4
            rad     = netsht + longw - 2*longout
         endif
      endif

* ----------------------------------------------------------------------
* calculate the vapor mass flux between the canopy and the surrounding 
* air mass
* ----------------------------------------------------------------------

      essnow = svp(tcanopy) 

* ----------------------------------------------------------------------
* added division by 10 to incorporate change in canopy resistance due
* to smoothing by intercepted snow
* ----------------------------------------------------------------------

      vaporflux = rhoair*(0.622/pgcm)*(vpair-essnow)/ra/10.0
      vaporflux = vaporflux/rhoh2o

      if (vpd.eq.0.and.vaporflux.lt.0.0) vaporflux = 0.0

* ----------------------------------------------------------------------
* calculate the latent heat flux 
* ----------------------------------------------------------------------

      ls = (677.-0.07*(tcanopy-t0c))*4.1868*1000.0
      lathsc = ls * vaporflux * rhoh2o

* ----------------------------------------------------------------------
* calculate the sensible heat flux
* in original c-vic there is no division of 10 after ra
* however, sometimes senhsc is too big, so the division is added.
* also, the lathsc and senhsc should not too high. so the maximum
* of lathsc and senhsc of 50 w/m2 is assigned.
* for completeness, we should use the energy balance mode for computing
* the canopy energy balance. june 2004/ji
* ----------------------------------------------------------------------
      senhsc = rhoair * cpair * (tair-tcanopy)/ra/10.0
  
* ----------------------------------------------------------------------
* calculate the advected energy
* ----------------------------------------------------------------------
      advectede = (ch_water*(tair-tcanopy)*rainf)/dtime

* ----------------------------------------------------------------------
* control the quality of computation
* ----------------------------------------------------------------------
      lathsc = max(min(lathsc, 50.0),-50.0)
      senhsc = max(min(senhsc, 50.0),-50.0)
      advectede = max(min(advectede, 50.0),-50.0)

* ----------------------------------------------------------------------
* calculate the amount of energy available for refreezing
* ----------------------------------------------------------------------
      refreezee = senhsc + lathsc + rad + advectede
      refreezee = refreezee * dtime             ! units in (j/m2)

* ----------------------------------------------------------------------
* if refreezee is positive it means energy is available to melt
* the intercepted snow in the canopy. if it is negative, it means that 
* intercepted water will be refrozen
* ----------------------------------------------------------------------

      vaporf = vaporflux * dtime   ! m/s -> m (+ air -> canopy)

      if (refreezee .gt. 0.0) then
         if (-(vaporf) .gt. cwater) then
            vaporf =  - cwater
            cwater = 0.0
         else
            cwater = cwater + vaporf
         end if
         
         potsnwmt = min((refreezee/lf/rhoh2o), csnow)
         
         maxwatint = liq_wc*csnow + maxint
         if ((cwater + potsnwmt) .le. maxwatint) then
            csnow  = csnow - potsnwmt
            cwater = cwater + potsnwmt
            potsnwmt = 0
         else
            exssnwmt = potsnwmt + cwater - maxwatint
            csnow = csnow + cwater - maxwatint
            cwater = maxwatint
            
            tmpintstr = 0.0
            releasmss = 0.0
            drip = 0.0
            
            if (snowtf.gt.0.0.and.inicsnow.le.min_is) then

* ----------------------------------------------------------------------
* water in excess of maxwatint has been generated.  if it is 
* snowing and there was little intercepted snow at the beginning 
* of the time step ( <= min_is), then allow the
* snow to melt as it is intercepted
* ----------------------------------------------------------------------

               drip  = exssnwmt
               csnow = csnow - exssnwmt
            else 
* ----------------------------------------------------------------------
* else, snowthroughfall = 0.0 or snowthroughfall > 0.0 and there is a 
* substantial amount of intercepted snow at the beginning of the time 
* step ( > min_is). snow melt may generate mass release.
* ----------------------------------------------------------------------

               tmpintstr = exssnwmt
            end if
            
            call massrelease(csnow,tmpintstr,releasmss,drip)

            raintf = raintf + drip
            snowtf = snowtf + releasmss
         end if

* ----------------------------------------------------------------------
* if intercepted snow has melted, add the water it held to drip
* ----------------------------------------------------------------------

         maxwatint = liq_wc*csnow + maxint
         if (cwater .gt. maxwatint) then
            drip = cwater - maxwatint
            cwater = maxwatint
            raintf = raintf + drip
         end if
  
      else                      ! else (refreezee <= 0.0) 
         
* ----------------------------------------------------------------------
* refreeze as much surface water as you can
* ----------------------------------------------------------------------

         if (refreezee .ge. (- cwater*lf))then
            csnow = csnow + abs(refreezee)/lf
            cwater = cwater - abs(refreezee)/lf
            refreezee = 0.0
         else 

* ----------------------------------------------------------------------
* all of the water in vegetation has been frozen.
* ----------------------------------------------------------------------
            csnow = csnow + cwater
            cwater = 0.0
         end if
         
         if (-(vaporf) .gt. csnow)then
            vaporf = - csnow
            csnow = 0.0
         else
            csnow = csnow + vaporf
         end if
      end if

      return
      end
