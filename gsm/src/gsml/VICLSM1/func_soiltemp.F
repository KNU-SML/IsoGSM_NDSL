	function func_soiltemp(t)

#include <vartyp.h>

************************************************************************
* source file:       func_soiltemp.f
* purpose:           soil temperature equation
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june & july 2003
* prgm history:      modified from 4.0.3 uw vic (soil_thermal_eqn.c)
*
* comments:
*      please see equation (8) in cherkauer and lettenmaier (jgr 1999)
*      for the soil thermal equation
************************************************************************

* ------------------- input variables ----------------------------------
	real t                  ! current node soil temperature (k)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
	real func_soiltemp      ! residure of thermal balance function
* ----------------------------------------------------------------------

* ------------------- common blocks ------------------------------------
#include <vic_soileb.h>
* ----------------------------------------------------------------------

* ------------------- local variables ----------------------------------
	real max_unfwat         ! function name
	real t0c		! ice/water mix temperature (k)
	parameter (t0c=273.15)
* ----------------------------------------------------------------------

	if(t.lt.t0c) then
	   ice = soilm - max_unfwat(t, soilmx, bubble, expt)
	else
	   ice = 0.0
	end if

	if(ice.lt.0.) ice = 0

	func_soiltemp = (t-t0c)*ee - aa*(tl-tu) - 
     &         bb*(tl+tu-2*t0c-gamma*fprime)-cc-dd*(ice-ice0)

#ifdef DBGVIC	
c	print *,'in func_soiltemp ',func_soiltemp,t,tl,tu,
c     &         gamma,fprime,ice,ice0,soilm,soilmx
c        print *,'   func_soiltemp ',(t-t0c)*ee,aa*(tl-tu),
c     &          bb*(tl+tu-2*t0c-gamma*fprime),cc,dd*(ice-ice0)
#endif
	return
	end
