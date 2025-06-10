      subroutine snow_melt(dtime,    vtype,   month, rhoair,
     &                        Ra,     wind,   vpair,    vpd,
     &                      Tair,    Tgrnd,  netsht,  longw,
     &                      pgcm,     frTF,    fsTF,  snowd,
     &                      frsn,      swq,     Tsf,    Tpk,
     &                   sfwater,  pkwater,   snowm,  latHs,
     &                     senHs)

#include <vartyp.h>

************************************************************************
* source file:       snow_melt.F
* purpose:           Calculate snow accumulation and melt using an energy 
*                    balance approach for a two layer snow model
* PRGMMR:            Ji Chen   
* ORG:               ECPC/CRD/SIO/UCSD 
* DATE:              JUNE, JULY, AUGUST & SEPTEMBER 2003
* PRGM HISTORY:      modified from 4.0.3 UW VIC (snow_melt.c)
************************************************************************

* ------------------- input variables ----------------------------------
 
* -- model basic parameter
      real dtime          ! time step (second)
      integer month       ! month of current time step

* -- atmosphere data

      real fsTF           ! Amount of snow (m/time step)
      real frTF           ! Rain through fall (m/time step)
      real Tair           ! Air temperature (K)
      real rhoair         ! air density  (kg/m3)
      real vpair          ! Actual vapor pressure of air (Pa)
      real vpd            ! Vapor pressure deficit (Pa)
      real Ra             ! aerodynamic resistance (s/m)
      real wind           ! Wind speed (m/s)
      real pgcm           ! Air pressure (Pa)

      real netsht         ! net shortwave radiation (w/m2)
      real longw          ! Longwave radiation (w/m2)

* -- vegetation parameters

      integer vtype       ! vegetation type

* -- land surface parameters

      real Tgrnd          ! ground surface temperature (K)
      real frsn           ! snow density (kg/m3)
      real snowd          ! snow depth (m)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real swq            ! Snow water equivalent at current pixel (m)
      real Tsf            ! Temperature of snow pack surface layer (K)
      real Tpk            ! Temperature of snow pack (K)
      real pkwater        ! Liquid water content of snow pack (m)
      real sfwater        ! Liquid water in snow surface (m)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real snowm          ! Amount of snowpack melt (m/time step)
      real latHs          ! Latent heat exchange at surface (W/m2)
      real senHs          ! Sensible heat exchange at surface (W/m2)
      real Refreeze       ! Refreeze energy (W/m2)
* ----------------------------------------------------------------------

* ------------------- common block -------------------------------------
#include <vic_snoweb.h>
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real Ice, Le, Ls
      real oldTsf           ! Snow Surface temp for previous time (K)
      real sfswq            ! Snow water equivalent on snow surface (m)
      real Pkswq            ! Snow pack snow water equivalent (m)
      real SfCC             ! Cold content of surface layer snow pack (J)
      real PkCC             ! Cold content of snow pack (J)
      real SnowTFCC         ! Cold content of new snowfall (J)
      real DeltPkswq        ! Change in snow water equivalent of pack (m)
      real DeltPkCC         ! Change in cold content of the pack
      real vaporf           ! vapor flux (m/time step) (+: to snow surface)

      real MaxLiqWat        ! Maximum liquid water content of pack (m)
      real Qnet             ! Net energy exchange at the surface (W/m2)
      real RefrozenW        ! Amount of refrozen water (m)
      real tmp              ! temporary variable

      character*80 ctrfct   ! determine function

      real func_snoweb      ! function name
      real root_brent       ! function name

      real rhoh2o           ! water density (kg/m^3)
      real Lf               ! Latent heat of freezing (J/kg) at 0C
      real t0c              ! ice/water mix temperature (K)
      real CH_ICE           ! Volumetric heat capacity (J/(m3*C)) of ice
      parameter (Lf=3.337e5, t0c=273.15, CH_ICE=2100.0e3,rhoh2o=1.e3)

      real SNOW_DT          ! snow surface temperature DT (K)
      real MAX_SS           ! MAX SURFACE SWE (first snow layer) (m)
      real LIQ_WC           ! LIQUID WATER CAPACITY (faction)
      parameter (SNOW_DT=1.0, MAX_SS=0.125, LIQ_WC=0.035)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* convert units ( K -> C)
* ----------------------------------------------------------------------

      Tair  = Tair - t0c             ! k -> C
      Tgrnd = Tgrnd - t0c            ! k -> C
      Tsf   = Tsf - t0c              ! k -> C
      Tpk   = Tpk - t0c              ! k -> C

      oldTsf= Tsf

* ----------------------------------------------------------------------
* Initialize snowpack variables
* ----------------------------------------------------------------------

      Ice = swq - pkwater - sfwater
      if(Ice .lt. 0) Ice = 0.0

* ----------------------------------------------------------------------
* Reconstruct snow pack  (MAX SURFACE SWE (first snow layer) = 0.125 m)
* ----------------------------------------------------------------------

      if (Ice .gt. MAX_SS) then
         Sfswq = MAX_SS
      else
         Sfswq = Ice
      end if

      Pkswq = Ice - Sfswq
  
* ----------------------------------------------------------------------
* Calculate cold contents
* ----------------------------------------------------------------------

      SfCC = CH_ICE * Sfswq * Tsf
      PkCC = CH_ICE * Pkswq * Tpk

      if(fsTF.gt.0) then            ! snowthrough fall
         if (Tair .gt. 0.0) then
            SnowTFCC = 0.0
         else
            SnowTFCC = CH_ICE * fsTF * Tair
         end if

* ----------------------------------------------------------------------
* Distribute fresh snowfall
* ----------------------------------------------------------------------

         if (fsTF .gt. (MAX_SS - Sfswq))then
            DeltPkswq = Sfswq + fsTF - MAX_SS
            
            if (DeltPkswq .gt. Sfswq) then
               DeltPkCC = SfCC + (fsTF - MAX_SS)/
     &                    fsTF * SnowTFCC
            else
               DeltPkCC = DeltPkswq/Sfswq * SfCC
            end if
            
            Sfswq = MAX_SS
            SfCC  = SfCC + SnowTFCC - DeltPkCC
            Pkswq = Pkswq + DeltPkswq
            PkCC  = PkCC + DeltPkCC
         else
            Sfswq = Sfswq + fsTF
            SfCC  = SfCC + SnowTFCC
         end if
      end if

      Tsf = SfCC/(CH_ICE * Sfswq)

      if (Pkswq .gt. 0.0) then   
         Tpk = PkCC/(CH_ICE * Pkswq)
      else
         Tpk = min(Tgrnd, min(0.0,Tair))
      end if

* ----------------------------------------------------------------------
* Adjust ice and snow surf_water
* ----------------------------------------------------------------------
      Ice = Ice + fsTF
      sfwater = sfwater + frTF

* ----------------------------------------------------------------------
* intialize vicsnoweb common block
* ----------------------------------------------------------------------
      Edtime  = dtime
      Emonth  = month
      EfrTF   = frTF
      ETair   = Tair
      Erhoair = rhoair
      Evpair  = vpair
      Evpd    = vpd
      ERa     = Ra
      Ewind   = wind
      Epgcm   = pgcm
      Enetsht = netsht
      Elongw  = longw
      Evtype  = vtype
      EoldTsf = oldTsf
      ETgrnd  = Tgrnd
      Efrsn   = frsn
      Esnowd  = snowd
      Esfswq  = sfswq
      Esfwater= sfwater

* ----------------------------------------------------------------------
* Calculate the snow surface energy balance for Tsf = 0.0 
* ----------------------------------------------------------------------

      Qnet = func_snoweb(0.0)

      Refreeze = Orefreeze
      latHs = OlatHs
      senHs = OsenHs

* ----------------------------------------------------------------------
* If Qnet is 0.0, then set the surface temperature to 0.0
* ----------------------------------------------------------------------

      if (abs(Qnet).lt.1.e-6)then
         Tsf = 0.0
         if (Refreeze .ge. 0.0) then    ! positive: frozen surface water
            RefrozenW = Refreeze/(Lf * rhoh2o)* dtime
            if (RefrozenW .gt. sfwater) then
               RefrozenW = sfwater
               Refreeze = RefrozenW * Lf * rhoh2o / dtime
            end if
            Sfswq   = Sfswq + RefrozenW
            Ice     = Ice + RefrozenW
            sfwater = sfwater - RefrozenW
            snowm   = 0.0
         else

* ----------------------------------------------------------------------
* Calculate snow melt
* ----------------------------------------------------------------------
            snowm  = abs(Refreeze)/(Lf *rhoh2o)*dtime      ! m/time step
         end if

* ----------------------------------------------------------------------
* Convert vapor mass flux to a depth per timestep and adjust surf_water
* ----------------------------------------------------------------------

         Le = 2.501e6
         vaporf = latHs/(Le * rhoh2o)

* ----------------------------------------------------------------------
* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990)
* ----------------------------------------------------------------------

         vaporf = vaporf * dtime   ! m/s -> m/time step (+: to snow surface)
        
         if (sfwater .lt. (-vaporf)) then
            vaporf  = -sfwater
            sfwater = 0.0
         else
            sfwater = sfwater + vaporf
         end if

* ----------------------------------------------------------------------
* If snowm < Ice, there was incomplete melting of the pack
* assume: melting snow from bottom to snow surface
* ----------------------------------------------------------------------

         if (snowm .lt. Ice) then
            if (snowm .le. Pkswq) then
               sfwater = sfwater + snowm
               Pkswq   = Pkswq - snowm
               Ice     = Ice   - snowm
            else 
               sfwater = sfwater + snowm + pkwater
               pkwater = 0.0
               Pkswq   = 0.0
               Ice     = Ice - snowm
               Sfswq   = Ice
            end if
    
* ----------------------------------------------------------------------
* Else, Snowm > Ice and there was complete melting of the pack
* ----------------------------------------------------------------------

         else
            snowm  = Ice
            sfwater= sfwater + Ice + pkwater
            Ice    = 0.0
            Sfswq  = 0.0
            Pkswq  = 0.0
            pkwater= 0.0
            Tsf    = min(0.0,Tair)
            Tpk    = min(0.0,Tair)
         end if

* ----------------------------------------------------------------------
* Else, SnowPackEB(T=0.0) <= 0.0 
* ----------------------------------------------------------------------

      else  

* ----------------------------------------------------------------------
* Calculate surface layer temperature using "Brent method"
* ----------------------------------------------------------------------

         ctrfct = 'func_snoweb'
         Tsf = root_brent(0.0, Tsf-SNOW_DT, ctrfct)

         Refreeze = Orefreeze
         latHs = OlatHs
         senHs = OsenHs

* ----------------------------------------------------------------------
* since we iterated, the surface layer is below freezing and no snowm
* ----------------------------------------------------------------------

         snowm  = 0.0
    
* ----------------------------------------------------------------------
* Since updated snow_temp < 0.0, all of the liquid water in the surface
* layer has been frozen   
* ----------------------------------------------------------------------

         SfSwq   = SfSwq + sfwater
         Ice     = Ice   + sfwater
         sfwater = 0.0
    
* ----------------------------------------------------------------------
* Convert vapor mass flux to a depth per timestep and adjust surf_water
* ----------------------------------------------------------------------

         Ls = (677. - 0.07 * Tsf) * 4.1868e3

* ----------------------------------------------------------------------
* Calculate latent heat flux
* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990)
* ----------------------------------------------------------------------
         vaporf = latHs/(Ls * rhoh2o)
         vaporf = vaporf * dtime      ! m/time step (+: to snow surface)

         if (Sfswq .lt. -vaporf) then ! in one time the maxi vaporf=Sfswq
            vaporf = -Sfswq
            Sfswq  = 0.0
            Ice    = Pkswq
         else
            Sfswq = Sfswq + vaporf
            Ice   = Ice + vaporf
         end if
      end if
  
* ----------------------------------------------------------------------
* Done with iteration etc, now Update the liquid water content of the
* surface layer  
* ----------------------------------------------------------------------

      MaxLiqWat = LIQ_WC * Sfswq

      if  (sfwater .gt. MaxLiqWat) then
         snowm = sfwater - MaxLiqWat
         sfwater  = MaxLiqWat
      else
         snowm = 0.0
      end if
  
* ----------------------------------------------------------------------
* Refreeze liquid water in the pack. variable 'Refreeze' is the
* heat released to the snow pack if all liquid water were refrozen.  
* if Refreeze < PackCC then all water IS refrozen  PackCC always <=0.0
* WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this
*  does not involve energy transported to the pixel. Instead heat from
*  the snow pack is used to refreeze water
* ----------------------------------------------------------------------
  
* ----------------------------------------------------------------------
* add surface layer outflow to pack liquid water
* ----------------------------------------------------------------------

      pkwater = pkwater + snowm

      Refreeze = pkwater * Lf * rhoh2o

* ----------------------------------------------------------------------
* calculate energy released to freeze
* ----------------------------------------------------------------------

      if (PkCC .lt. -Refreeze) then      ! cold content not fully depleted
         Pkswq = Pkswq + pkwater         ! refreeze all water and update
         Ice   = Ice + pkwater
         pkwater = 0.0
         if (Pkswq .gt. 0.0) then
c            PkCC = Pkswq * CH_ICE * Tpk + Refreeze ! original vic equation
            PkCC = PkCC + Refreeze
            Tpk  = PkCC / (CH_ICE * Pkswq)
         else 
            Tpk = 0.0
         end if
      else  

* ----------------------------------------------------------------------
* cold content has been either exactly satisfied or exceeded. If
* PackCC = refreeze then pack is ripe and all pack water is
* refrozen, else if energy released in refreezing exceeds PackCC 
* then exactly the right amount of water is refrozen to satify PackCC.
* The refrozen water is added to PackSwq and Ice
* ----------------------------------------------------------------------
         Tpk      = 0.0
         DeltPkswq = -PkCC/(Lf * rhoh2o)  ! the amount of pkwater -> Pkswq

         pkwater = pkwater - DeltPkSwq
         Pkswq = Pkswq + DeltPkSwq
         Ice = Ice + DeltPkSwq 
      end if
  
* ----------------------------------------------------------------------
* Update the liquid water content of the pack
* ----------------------------------------------------------------------

      MaxLiqWat = LIQ_WC * Pkswq
      if (pkwater .gt. MaxLiqWat) then
         snowm = pkwater - MaxLiqWat
         pkwater = MaxLiqWat
      else
         snowm = 0.0
      end if

* ----------------------------------------------------------------------
* Update snow properties
* ----------------------------------------------------------------------
      Ice  = Pkswq + Sfswq

      if (Ice .gt. MAX_SS) then
         SfCC   = CH_ICE * Tsf * Sfswq
         PkCC   = CH_ICE * Tpk * Pkswq
         if (Sfswq .gt. MAX_SS) then
            tmp = SfCC * (Sfswq-MAX_SS)/Sfswq
            PkCC  = PkCC + tmp
            SfCC  = SfCC - tmp
            Pkswq = Pkswq + Sfswq - MAX_SS
            Sfswq = MAX_SS
         else
            tmp = PkCC * (MAX_SS - Sfswq) / Pkswq
            PkCC  = PkCC - tmp
            SfCC  = SfCC + tmp
            Pkswq = Pkswq + Sfswq - MAX_SS
            Sfswq = MAX_SS
         end if

         Tpk = PkCC / (CH_ICE * Pkswq)
         Tsf = SfCC / (CH_ICE * Sfswq)
      else 
         Pkswq = 0.0
         Tpk = 0.0
      end if

      swq = Ice + pkwater + sfwater

* ----------------------------------------------------------------------
* convert units ( c -> K)
* ----------------------------------------------------------------------

      Tair  = Tair + t0c             ! c -> k
      Tgrnd = Tgrnd + t0c            ! c -> k
      Tsf   = Tsf + t0c              ! c -> k
      Tpk   = Tpk + t0c              ! c -> k

      return
      end
