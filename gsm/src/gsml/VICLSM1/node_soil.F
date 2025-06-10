      subroutine node_soil(msl,    nsl,    dphnd,    dph,
     &                   smxnd,  expnd,    bubnd,    sld,
     &                     bkd,    qrt,     fsmc,     tn,
     &                    smnd,  icend,    kapnd,   csnd)

#include <vartyp.h>

************************************************************************
* source file:       node_soil.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic c version
*             (distribute_node_moisture_properties in soil_conduction.c)
*
*  abstract: this subroutine determines the moisture and ice contents
*        of each soil thermal node based on the current node temperature
*        and layer moisture content. thermal conductivity and volumetric
*        heat capacity are then estimated for each node based on the
*        division of moisture contents. soil thermal conductivity
*        calculated using johansens method.
*
*  reference: farouki, o.t., "thermal properties of soils" 1986
*  chapter 7: methods for calculating the thermal conductivity of soils
*     h.b.h. - refers to the handbook of hydrology.
*
*  program history: node_soil.f is originally a part of 
*                   uw_vic (distribute_node_moisture_properties in 
*                           soil_conduction.c)
************************************************************************

* ------------------- input variables ----------------------------------
      integer msl        ! number of soil moisture layers
      integer nsl        ! number of soil thermal nodes
      real dphnd(nsl)    ! thermal node thickness (m)
      real dph(msl)      ! soil layer thickness (m)
      real smxnd(nsl)    ! thermal node maximum moisture content (mm/mm)
      real expnd(nsl)    ! thermal node exponential (n/a)
      real bubnd(nsl)    ! thermal node bubbling pressure (cm)
      real sld(msl)      ! soil density (kg m-3)
      real bkd(msl)      ! soil layer bulk density (kg m-3)
      real qrt(msl)      ! soil qrt content (fraction)
      real fsmc(msl)     ! soil layer moisture (mm)
      real tn(nsl)       ! thermal node temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real smnd(nsl)     ! thermal node moisture content (mm/mm)
      real icend(nsl)    ! thermal node ice content (mm/mm)
      real kapnd(nsl)    ! thermal node thermal conductivity (w m-1 k-1)
      real csnd(nsl)     ! thermal node heat capacity (j m-3 k-1)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real lsum, zsum

      real max_unfwat, soil_condvty, heat_capcty ! function names

      integer lidx, nidx
      logical past_bottom

      real t0c
      parameter (t0c = 273.15)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* initialize the computation
* ----------------------------------------------------------------------

      past_bottom = .false.
      lidx = 1
      lsum = 0.
      zsum = 0.

      do nidx = 1, nsl
         zsum = zsum + dphnd(nidx)
         if(zsum.gt.lsum.and..not.past_bottom) then
            lsum = lsum + dph(lidx)
            lidx = lidx + 1
            if( lidx .ge. msl ) then
               past_bottom = .true.
               lidx = msl
            end if

            do while(zsum.gt.lsum)
               lsum = lsum + dph(lidx)
               lidx = lidx + 1
               if( lidx .ge. msl ) then
                  past_bottom = .true.
                  lidx = msl
               end if
            end do
         end if

         if(zsum.eq.lsum.and.nidx.ne.1  ! node on layer boundary
     &        .and.lidx.ne.1)then
            smnd(nidx) = (fsmc(lidx-1)/dph(lidx-1) 
     &	                + fsmc(lidx)/dph(lidx))/1000.0/2.
         else                            ! node completely in layer
            smnd(nidx) = fsmc(lidx)/dph(lidx)/1000.0
         end if

         if(tn(nidx).lt.t0c) then        ! compute moisture/ice contents
            icend(nidx) = smnd(nidx) - 
     &        max_unfwat(tn(nidx),smxnd(nidx),
     &                   bubnd(nidx),expnd(nidx))

            if(icend(nidx).lt.0) icend(nidx)=0

* ----------------------------------------------------------------------
* compute thermal conductivity
* ----------------------------------------------------------------------
            kapnd(nidx) = soil_condvty(smnd(nidx),
     &                smnd(nidx)-icend(nidx),
     &                sld(lidx),bkd(lidx),qrt(lidx))
         else                          !compute moisture and ice contents
            icend(nidx)   = 0
* ----------------------------------------------------------------------
* compute thermal conductivity
* ----------------------------------------------------------------------
            kapnd(nidx)= soil_condvty(smnd(nidx), 
     &                smnd(nidx), sld(lidx), bkd(lidx), qrt(lidx))
         end if

* ----------------------------------------------------------------------
* compute volumetric heat capacity 
* ----------------------------------------------------------------------

         csnd(nidx) = heat_capcty(bkd(lidx)/sld(lidx),
     &                smnd(nidx)-icend(nidx),icend(nidx))
      end do

      return
      end 
