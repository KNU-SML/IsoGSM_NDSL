      subroutine rainorsnow(msl,   nsl,  msub,    md,
     2                     nveg,   prc,  tgcm,   fmu,
     3                    lstsn, rfall, sfall,   cwt,
     4                      csn,   swq,   rsn,   tsf,
     5                      tpk,   sfw,   pkw,   smc,
     6                      sic,   tnd)

#include <vartyp.h>

************************************************************************
* source file:       rainorsnow.F
* PRGMMR:            Ji Chen   
* ORG:               ECPC/CRD/SIO/UCSD 
* DATE:              JUNE, JULY, AUGUST & SEPTEMBER 2003
* program history:   modified from uw_vic (calc_rainonly.c,dist_prec.c,...)
*
* subroutine documentation block
*  
* abstract: decides rainfall or snowfall, and redistribute hydrologic
*           properties in the wet and dry districts according to 
*           rainfall intensity
*
* comments: if snow falling, mu = 1.0. if no snow falling, mu varies
*           according to the rainfall intensity, and the related 
*           hydrologic variables, such as csnow, cwater, raintf, snowtf, 
*           and other soil variables are recalcuated. if no precipitation,
*           mu doesnt change and its value is decided by the last 
*           precipitation event. 
************************************************************************

* ------------------- input variables ----------------------------------
      integer msl           ! number of soil layer
      integer nsl           ! number of soil node
      integer msub          ! maximum number of subgrids
      integer md            ! number of precipitation district
      integer nveg          ! number of subgrids

      real prc              ! precipitation (rain+snow) (m/time step)
      real tgcm             ! air temperature (k) 
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fmu              ! fraction of wet district
      real cwt(md,msub)     ! canopy intercepted water 1:wet, 2:dry (m)
      real csn(md,msub)     ! canopy intercepted snow  (m)
      real swq(md,msub)     ! snow water equivalent (m)
      real rsn(md,msub)     ! snow density (kg/m3)
      real tsf(md,msub)     ! snow surface temperature (k)
      real tpk(md,msub)     ! snow pack temperature (k)
      real sfw(md,msub)     ! snow surface water (m)
      real pkw(md,msub)     ! snow pack water (m)
      real smc(msl,md,msub) ! soil moisture (mm)
      real sic(msl,md,msub) ! soil ice (mm)
      real tnd(nsl,md,msub) ! soil temperature (k)

      integer lstsn         ! time steps since last snowfall (time step)
                            ! lstsn = -1 (no snow), 0 (new snow)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real rfall            ! rainfall (m/time step)
      real sfall            ! snowfall (m/time step)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real ta               ! air temperature (c)
      real oldmu, oldwu, wu

      integer nv, n, m, idi ! loop indexes
      logical logicswq

      real tmaxsnow         ! maxi temperature at which snow can fall (c)
      real tminrain         ! mini temperature at which rain can fall (c)
      real t0c              ! ice/water mix temperature (k)
      real prec_expt        ! exponential that controls the fraction of a
                            ! grid cell that receives rain during a storm
                            ! of given intensity

      parameter (tmaxsnow=0.5, tminrain=-0.5)
      parameter (t0c=273.15, prec_expt = 0.6)
* ----------------------------------------------------------------------

      if(prc.gt.0.0) then

         oldmu = fmu
         oldwu = 1.0 - oldmu

         ta = tgcm - t0c

         if(ta.lt.tmaxsnow.and.ta.gt.tminrain)then
            rfall = (ta-tminrain)/(tmaxsnow-tminrain)*prc
         else if(ta .gt.tmaxsnow) then
            rfall = prc
         end if

         sfall = prc - rfall
         
* ----------------------------------------------------------------------
* the fractional coverage of precipitation over an area or grid cell, 
*     mu, is estimated using the equation from fan et. al. (wrr 1996)
*     if snowfalling, mu = 1.0
* ----------------------------------------------------------------------

         if(sfall.gt.0.0) then
            lstsn = 0
            fmu = 1.0
         else
            if(md.eq.1) then
               fmu = 1.0
            else
               fmu = 1.0 - exp(-prec_expt*rfall) 
            end if
         end if
         
* ----------------------------------------------------------------------
* reorganize variables in the wet and dry districts for continous run
* ----------------------------------------------------------------------

         if(md.ne.1) then      ! include precipitation distribution

            if(fmu.ne.oldmu) then
               wu = 1.0 - fmu

               do nv = 1, nveg
                  if(fmu.gt.oldmu) then             ! part of dry -> wet
                     if(csn(1,nv).gt.0.or.csn(md,nv).gt.0) then
                        csn(1,nv) = (csn(1,nv)*oldmu+
     &                              csn(md, nv)*(fmu-oldmu))/fmu
                     end if
                     if(cwt(1,nv).gt.0.or.cwt(md,nv).gt.0) then
                        cwt(1,nv)= (cwt(1,nv)*oldmu+
     &                             cwt(md, nv)*(fmu-oldmu))/fmu
                     end if
                  else                              ! part of wet -> dry
                     if(csn(1,nv).gt.0.or.csn(md,nv).gt.0) then
                        csn(md,nv) = (csn(md,nv)*oldwu+
     &                              csn(1, nv)*(wu-oldwu))/wu
                     end if
                     if(cwt(1,nv).gt.0.or.cwt(md,nv).gt.0) then
                        cwt(md,nv)= (cwt(md,nv)*oldwu+
     &                             cwt(1, nv)*(wu-oldwu))/wu
                     end if
                  end if

                  if(swq(1,nv).gt.0.or.swq(md,nv).gt.0) then
                     if(fmu.gt.oldmu) then
                        swq(1,nv) = (swq(1,nv)*oldmu+
     &                              swq(md,nv)*(fmu-oldmu))/fmu
                        rsn(1,nv) = (rsn(1,nv)*oldmu+
     &                              rsn(md,nv)*(fmu-oldmu))/fmu
                        tsf(1,nv) = (tsf(1,nv)*oldmu+
     &                              tsf(md,nv)*(fmu-oldmu))/fmu
                        tpk(1,nv) = (tpk(1,nv)*oldmu+
     &                              tpk(md,nv)*(fmu-oldmu))/fmu
                        sfw(1,nv) = (sfw(1,nv)*oldmu+
     &                              sfw(md,nv)*(fmu-oldmu))/fmu
                        pkw(1,nv) = (pkw(1,nv)*oldmu+
     &                              pkw(md,nv)*(fmu-oldmu))/fmu
                     else
                        swq(md,nv) = (swq(md,nv)*oldwu+
     &                              swq(1,nv)*(wu-oldwu))/wu
                        rsn(md,nv) = (rsn(md,nv)*oldwu+
     &                              rsn(1,nv)*(wu-oldwu))/wu
                        tsf(md,nv) = (tsf(md,nv)*oldwu+
     &                              tsf(1,nv)*(wu-oldwu))/wu
                        tpk(md,nv) = (tpk(md,nv)*oldwu+
     &                              tpk(1,nv)*(wu-oldwu))/wu
                        sfw(md,nv) = (sfw(md,nv)*oldwu+
     &                              sfw(1,nv)*(wu-oldwu))/wu
                        pkw(md,nv) = (pkw(md,nv)*oldwu+
     &                              pkw(1,nv)*(wu-oldwu))/wu
                     end if
                  end if

                  do m = 1, msl
                     if(fmu.gt.oldmu) then             ! part dry -> wet
                        smc(m,1,nv) = (smc(m,1,nv)*oldmu+
     &                                smc(m,md,nv)*(fmu-oldmu))/fmu
                        sic(m,1,nv) = (sic(m,1,nv)*oldmu+
     &                                sic(m,md,nv)*(fmu-oldmu))/fmu
                     else                              ! part wet -> dry
                        smc(m,md,nv) = (smc(m,md,nv)*oldwu+
     &                                smc(m,1,nv)*(wu-oldwu))/wu
                        sic(m,md,nv) = (sic(m,md,nv)*oldwu+
     &                                sic(m,md,nv)*(wu-oldwu))/wu
                     end if
                  end do

                  do n = 1, nsl
                     if(fmu.gt.oldmu) then
                        tnd(n,1,nv) = (tnd(n,1,nv)*oldmu + 
     &                                tnd(n,md,nv)*(fmu-oldmu))/fmu
                     else
                        tnd(n,md,nv) = (tnd(n,md,nv)*oldwu + 
     &                                tnd(n,1,nv)*(wu-oldwu))/wu
                     end if
                  end do
               end do           ! loop of nveg
            end if
         end if             ! end of precip distribution
      end if           ! end of prc>0
         
      if(sfall.le.0.0) then
         logicswq = .false.
         do nv = 1, nveg         !whether or not snow exists on ground
            if(swq(1,nv).gt.0.or.swq(md,nv).gt.0) then
               logicswq = .true.
            end if
         end do

         if(logicswq) then
            if(lstsn.ge.0) then
               lstsn = lstsn + 1
            else
               lstsn = 0
            end if
         else
            lstsn = -1
         end if
      end if

      return
      end 
