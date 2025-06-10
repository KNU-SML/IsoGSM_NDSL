      function soil_condvty(moist, wu, rhosoil, rhobulk,
     &                          quartz)

#include <vartyp.h>

************************************************************************
* source file:       soil_condvty.f
* purpose:           soil thermal conductivity using johansens method
* prgmmr:            ji chen
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic 
*                    (soil_conductivity in soil_conduction.c)
*
*  reference: farouki, o.t., "thermal properties of soils" 1986
*     chapter 7: methods for calculating the thermal conductivity of
*     soils.     h.b.h. - refers to the handbook of hydrology.
*
*  porosity = n = porosity
*  ratio = sr = fractionaldegree of saturation
*  all k values are conductivity in w/mk
*  wu is the fractional volume of unfrozen water
*
************************************************************************

* ------------------- input variables ----------------------------------
      real moist         ! total moisture content (mm/mm)
      real wu            ! liquid water content (mm/mm)
      real rhosoil       ! soil density (kg m-3)
      real rhobulk       ! soil bulk density (kg m-3)
      real quartz        ! soil quartz content (fraction)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real soil_condvty  ! soil thermal conductivity (w m-1 k-1)
*-----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real ke
      real ki            ! thermal conductivity of ice (w/mk)
      real kw            ! thermal conductivity of water (w/mk)
      real ksat
      real ks	         ! thermal conductivity of solid (w/mk)
			 ! function of quartz content
      real kdry
      real sr	         ! fractional degree of saturation
      real k
      real porosity
      parameter (ki=2.2, kw = 0.57)
* ----------------------------------------------------------------------

      kdry = (0.135*rhobulk+64.7)/(rhosoil - 0.947*rhobulk)

      if(moist.gt.0.) then
         porosity = 1.0 - rhobulk / rhosoil
         sr = moist/porosity
         ks = (7.7**quartz) * (2.2**(1.0-quartz))

         if(wu.eq.moist) then             !soil unfrozen
            ksat = (ks**(1.0-porosity)) * (kw**porosity)
            ke = 0.7 * log10(sr) + 1.0
         else                             !soil frozen
            ksat = (ks**(1.0-porosity))*(ki**(porosity-wu))*(kw**wu)
            ke = sr
         end if

         k = (ksat-kdry)*ke+kdry
         if(k.lt.kdry) k=kdry
      else
         k = kdry
      end if

      soil_condvty = k

      return
      end
