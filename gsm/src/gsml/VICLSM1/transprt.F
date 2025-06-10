      subroutine transprt(msl, dtime,  fsmc,    fsic,
     &                    wcr,   wpw,  root,       f,
     &                   frac,   rad,   vpd,      ra,
     &                     rs,  rarc, fflai,  netsht,
     &                   tair,  pgcm,   rgl,   mfsmc,
     &                 layevp)

#include <vartyp.h>

************************************************************************
* source file:       transprt.f
* purpose:           computes the traspiration
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (canopy_evap.c)
*
* comments:
*    please check equations (5), (7), (8), (9), (10), and (11) in liang, 
*    lettenmaier, wood, and burges (wrr 1994) for computing evap.
************************************************************************

* ------------------- input variables ----------------------------------

      integer msl      ! number of soil layer
      real dtime       ! time step (second)
      real fsic(msl)   ! soil ice (mm)
      real wcr(msl)    ! ~70% field capacity (mm)
      real wpw(msl)    ! wilting point soil moisture (mm)
      real root(msl)   ! root content (fraction)
      real f           ! fraction of time step for no canopy water
      real frac        ! fraction of wet leaf
      real rad         ! net radiation (w/m2)
      real vpd         ! vapor pressure deficit (pa)
      real ra          ! aerodynamical resistance (s/m)
      real rs          ! minimum stomatal resistance (s/m)
      real rarc        ! architectural resistance (s/m)
      real fflai       ! leaf area index
      real netsht      ! net shortwave radiation (w/m2)
      real tair        ! air temperature (k)
      real pgcm        ! pressure (pa)
      real rgl         ! value of solar radiation below which there 
                       ! will be no transpiration (ranges from ~30 w/m^2
                       ! for trees to ~100 w/m^2 for crops)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fsmc(msl)   ! soil moisture (liquid + ice) (mm)
      real mfsmc(msl)  ! soil moisture (liquid + ice) (mm)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real layevp(msl) ! transpiration from each soil layer (m/time step)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      integer k        ! loop index
      real gsm_inv     ! soil moisture stress factor
      real evap        ! tmp holding for evap total
      real moist1      ! tmp holding of moisture top
      real moist2      ! tmp holding of moisture bottom
      real wcr1        ! tmp holding of critical water top
      real wcr2        ! tmp holding of critical water bottom
      real root1       ! tmp holding of root top
      real root2       ! tmp holding of root bottom
      real rootsum     ! proportion of roots in moist>wcr zones
      real spare_evap  ! evap for 2nd distribution
      real avlsm(msl)  ! moisture available for trans (mm)

      real vicpenman   ! function name
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* computes evapotranspiration for unfrozen soils. allows multiple layers
* compute moisture content in combined upper layers
* ----------------------------------------------------------------------

      moist1 = 0.0
      wcr1 = 0.0
      root1= 0.0

      do k = 1, msl-1
         if(root(k) .gt. 0.) then
            avlsm(k) = fsmc(k)-fsic(k)
            moist1= moist1 + avlsm(k)
            wcr1  = wcr1 + wcr(k)
            root1 = root1 + root(k)
         else
            avlsm(k) = 0.0
         end if
      end do

* ----------------------------------------------------------------------
* compute moisture content in lowest layer
* ----------------------------------------------------------------------

      k = msl
      moist2   = fsmc(k)-fsic(k)
      avlsm(k) = moist2
      wcr2  = wcr(k)
      root2 = root(k)

* ----------------------------------------------------------------------
* please refer the content (p14417-14418) of liang et al. (jgr 1994)
* for the following two cases of transpiration 
*
* case 1: moisture in both layers exceeds wcr, or moisture in
*     layer with more than half of the roots exceeds wcr.
*     potential evapotranspiration not hindered by soil dryness.  if
*     layer with less than half the roots is dryer than wcr, extra
*     evaporation is taken from the wetter layer. otherwise layers
*     contribute to evapotransipration based on root fraction.
* ----------------------------------------------------------------------

      if( (moist1.ge.wcr1.and.moist2.ge.wcr2.and.wcr1.gt.0.).or.
     &    (moist1.ge.wcr1.and.root1.ge.0.5) .or.
     &    (moist2.ge.wcr2.and.root2.ge.0.5) ) then

         gsm_inv=1.0

* -- evap (m/s)

         evap = vicpenman(rad,    vpd,      ra,      rs, 
     &                   rarc,  fflai, gsm_inv,    tair, 
     &                 netsht,   pgcm,     rgl)

         evap = evap*(1.0-f*frac)*dtime     ! m/s -> m/time step

* ----------------------------------------------------------------------
* divide up evap based on root distribution
* note the indexing of the roots 
* ----------------------------------------------------------------------

         rootsum = 1.0
         spare_evap = 0.0

         do k = 1, msl
            if(avlsm(k).ge.wcr(k))then
               layevp(k) = evap*root(k)
            else
               if (avlsm(k) .ge. wpw(k)) then 
                  gsm_inv = (avlsm(k) - wpw(k)) /
     &                 (wcr(k) - wpw(k))
               else 
                  gsm_inv = 0.0
               end if

               layevp(k)  = evap*gsm_inv*root(k)
               rootsum    = rootsum - root(k)
               spare_evap = evap*root(k)*(1.0-gsm_inv)
            end if
         end do

* ----------------------------------------------------------------------
* assign excess evaporation to wetter layer
* ----------------------------------------------------------------------

         if(spare_evap.gt.0.0) then
            do k = 1, msl
               if(avlsm(k).ge.wcr(k)) then
                  layevp(k) = layevp(k)+root(k)*spare_evap/rootsum
               end if
            end do
         end if
        
#ifdef DBGVIC
         print 77,evap,layevp
 77      format('case 1 tra evap',e10.4,' layevp',3(e10.4,1x))
#endif
* ----------------------------------------------------------------------
* case 2: independent evapotranspirations
*    evapotranspiration is restricted by low soil moisture. evaporation
*    is computed independantly from each soil layer.
* ----------------------------------------------------------------------

      else 
         do k = 1, msl

* ----------------------------------------------------------------------
* set evaporation restriction factor 
* ----------------------------------------------------------------------
            if(avlsm(k) .ge. wcr(k)) then
               gsm_inv = 1.0
            else if(avlsm(k) .ge. wpw(k)) then
               gsm_inv = (avlsm(k) - wpw(k)) /
     &                   (wcr(k) - wpw(k))
            else 
               gsm_inv = 0.0
            end if

            if(gsm_inv .gt. 0.0 .and. root(k).gt.0.0) then

* ----------------------------------------------------------------------
* compute potential evapotranspiration
* ----------------------------------------------------------------------

               evap = vicpenman(rad,   vpd,      ra,    rs,
     &                         rarc, fflai, gsm_inv,  tair,
     &                       netsht,  pgcm,     rgl)

               layevp(k) = evap*root(k)*(1.0-f*frac)*dtime
#ifdef DBGVIC
               print 66,evap,k,layevp(k)
 66            format('case 2 evap',e10.4,' k',i4,' layevp',e10.4)
#endif
            else
               layevp(k) = 0.0
            end if
         end do
      end if

* ----------------------------------------------------------------------
* check that transpiration does not cause soil moisture to 
*       fall below wilting point.
* ----------------------------------------------------------------------

#ifdef DBGVIC
      print *,' bef fsmc-wpw layevp',layevp
#endif

      do k = 1, msl
         if(layevp(k)*1000.0 .gt. (fsmc(k)-wpw(k))) then
            layevp(k) = (fsmc(k) - wpw(k))/1000.0
         endif
         if (layevp(k) .lt. 0.0 ) then
            layevp(k) = 0.0
         end if
      end do

#ifdef DBGVIC
      print *,' aft fsmc-wpw layevp',layevp
#endif

      do k = 1, msl
         mfsmc(k) = fsmc(k) - layevp(k)*1000.0
      end do

      return
      end
