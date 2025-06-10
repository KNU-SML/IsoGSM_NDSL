      function canopy_evap(msl,  vtype,   dtime,  fflai,
     &                    fcwt,  rfall,    fsmc,   fsic,
     &                     wcr,    wpw,     frt,    rad,
     &                  rhoair,  vpair,     vpd,     ra,
     &                  netsht,   tair,    pgcm,  mfcwt,
     &                   mfrtf,  mfsmc)

#include <vartyp.h>

************************************************************************
* source file:       canopy_evap.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (canopy_evap.c)
*
* abstract: this function computes the evaporation, traspiration and 
*           throughfall of the vegetation types for multi-layered model.
*
* comments:
*    please check equations (1), (5), (7), (8), (9), (10), and (11) 
*    in liang, lettenmaier, wood, and burges (wrr 1994) for computing et
************************************************************************

* ------------------- input variables ----------------------------------
      integer msl       ! number of soil layer
      integer vtype     ! vegetation type
      real dtime        ! time step (second)
      real fflai        ! leaf area index
      real fsic(msl)    ! soil ice (mm)
      real wcr(msl)     ! ~70% field capacity (mm)
      real wpw(msl)     ! wilting point (mm)
      real frt(msl)     ! root content (fraction)
      real rad          ! net radiation (w/m2)
      real rhoair       ! air density (kg/m^3)
      real vpair        ! actual vapor pressure of air (pa)
      real vpd          ! vapor pressure deficit (pa)
      real ra           ! aerodynamical resistance (s/m)
      real netsht       ! net shortwave radiation (w/m2)
      real tair         ! air temperature (k)
      real pgcm         ! pressure (pa) 
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real rfall        ! precipitation during time step (m/time step)
      real fcwt         ! canopy water (m)
      real fsmc(msl)    ! soil moisture (liquid + ice) (mm)
      real mfrtf        ! rain throughfall (m/time step)
      real mfcwt        ! canopy water (m)
      real mfsmc(msl)   ! soil moisture (liquid + ice) (mm)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real canopy_evap  ! et (m/s)
* ----------------------------------------------------------------------

#include <vic_veglib.h>

* ------------------- local variables ----------------------------------
      integer k 

      real wdmax        ! maximum canopy water hold capacity (m)
      real evap         ! evaporation (m/s) (from vicpenman)
      real tmp_wdew     ! temporary trap canopy water (m)
      real frac         ! wet leaf fraction (fraction)
      real canopyevap   ! canopy evaporation from wet leaf (m/time step)

      real rs           ! minimum stomatal resistance (s/m)
      real rarc         ! architectural resistance (s/m)
      real rgl          ! value of solar radiation below which there 
                        ! will be no transpiration (ranges from ~30 w/m^2
                        ! for trees to ~100 w/m^2 for crops)

      real layevp(msl)  ! transpiration from each soillayer (m/time step)
      real vicpenman    ! function name
      real svp          ! function name

      real tmp_evap     ! temporary total evaporation (m/time step)
      real bas, f

      real rhoh2o       ! water density (kg/m^3)
      real lai_wf       ! leaf water factor for interception storage (m)
      parameter (rhoh2o=1.e3,lai_wf=0.0002)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* compute evaporation from canopy intercepted water
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* definitions:
*    wdmax - max dew and rain holding capacity
*    fcwt - dew and rain trapped on vegetation
* ----------------------------------------------------------------------

      rs  = veg_rmin(vtype)
      rarc= veg_rarc(vtype)
      rgl = veg_rgl(vtype)
      mfrtf = 0.0

      do k = 1, msl
         mfsmc(k) = fsmc(k)
      enddo

      tmp_wdew = fcwt
      wdmax = lai_wf * fflai

      if (tmp_wdew.gt.wdmax) then
         mfrtf = tmp_wdew - wdmax
         tmp_wdew = wdmax
      end if

#ifdef DBGVIC
      print *,'in canopy_evap',fcwt,wdmax,tmp_wdew
#endif

* ----------------------------------------------------------------------
* in the original c-vic (uw), the negative rad is also used in penman
* method, which will introduce errors in computing energy balance.
* therefore, we change the negative rad computation. if rad is negative,
* we use aerodynamic method to compute evap and assume canopy 
* temperature equals air temperature.
* ----------------------------------------------------------------------
* -- evap (m/s)

      if(rad.gt.0.0) then
         evap = vicpenman(rad,    vpd,    ra,    rs,
     &                   rarc,  fflai,   1.0,  tair,
     &                 netsht,   pgcm,   rgl)
      else
         evap = rhoair*(0.622/pgcm)*(svp(tair)-vpair)/ra
         evap = evap / rhoh2o
         if(vpd.gt.0.0.and.evap.lt.0.0) evap = 0.0  ! no dew
         if (vpd.eq.0.0.and.evap.gt.0.0) evap=0.0
         if (evap.gt.0.) evap = 0.0       ! negative rad, no evap
      endif

      if(evap.ge.0.0)then
         if(tmp_wdew.gt.0) then
            bas  = tmp_wdew / wdmax
            frac = bas**(2.0/3.0) ! wet leaf fraction
         else
            frac = 0.0
         endif
         canopyevap = frac * evap * dtime            
      else
         frac = 1.0
         canopyevap = evap*dtime     ! dew
      endif
         
* ----------------------------------------------------------------------
* evap can not exceed current storage
* ----------------------------------------------------------------------

      if (canopyevap .gt. 0.0) then
         f = min(1.0,(tmp_wdew / canopyevap))
      else
         f = 1.0
      endif
      canopyevap = canopyevap*f

      tmp_wdew = tmp_wdew + rfall - canopyevap
      if (tmp_wdew .gt. wdmax)then
         mfrtf  = mfrtf + tmp_wdew - wdmax
         tmp_wdew = wdmax
      endif

* ----------------------------------------------------------------------
* compute transpiration from vegetation
* ----------------------------------------------------------------------
      tmp_evap = canopyevap

* if evap .lt. 0, dewing occurs and the whole leaf wet and no transprt

#ifdef DBGVIC
      print 77,canopyevap,evap,rad,tmp_wdew,mfrtf
 77   format('canpevp=',e12.4,4(1x,e12.4))
#endif

* netsht radiation > 0 for activing transpiration processes (6/2004)
* ra/10.0 (6/2004)

      if(evap.ge.0.and.rad.ge.0.0.and.netsht.gt.0.0) then
         call transprt(msl,  dtime,   fsmc,    fsic,
     &                 wcr,    wpw,    frt,       f,
     &                frac,    rad,    vpd,      ra,
     &                  rs,   rarc,  fflai,  netsht,
     &                tair,   pgcm,    rgl,   mfsmc,
     &              layevp)

         do k = 1, msl
            tmp_evap = tmp_evap + layevp(k)
         end do

#ifdef DBGVIC
         print *,'transprt ',msl,  dtime,   fsmc,    fsic,
     &                 wcr,    wpw,    frt,       f,
     &                frac,    rad,    vpd,      ra,
     &                  rs,   rarc,  fflai,  netsht,
     &                tair,   pgcm,    rgl,   mfsmc,
     &              layevp
#endif

      endif

      mfcwt = tmp_wdew

      canopy_evap = tmp_evap/dtime     ! m/time step -> m/s

      return
      end 
