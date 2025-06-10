      subroutine aerodynamic(vtype,   month,   slrough, snrough,
     &                        zref,      wind,       u,      ra)

#include <vartyp.h>

************************************************************************
* source file:       aerodynamic.f
* prgmmr:            ji chen
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august, & sept. 2003
* prgm history:      modified from 4.0.3 uw vic (calcaerodynamic.c)
*
* abstract: calculate the aerodynamic resistance for each vegetation 
*     layer, and the wind 2m above the layer boundary.  in case of 
*     an overstory, also calculate the wind in the overstory.
*     the values are normalized based on a reference height wind 
*     speed, uref, of 1 m/s.  to get wind speeds and aerodynamic 
*     resistances for other values of uref, you need to multiply 
*     the here calculated wind speeds by uref and divide the 
*     here calculated aerodynamic resistances by uref
*   
* comments: 
*     please study the paper by wigmosta, vail and lettenmaier (wrr 1994)
*     for understanding the equations used to compute the ra
************************************************************************

* ------------------- input variables ----------------------------------
      integer vtype       ! vegetation type
      integer month       ! month of current time step
      real slrough        ! soil roughness (m)
      real snrough        ! snow roughness (m)
      real zref           ! reference height for windspeed (m)
      real wind           ! wind speed at zref (m/s)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
* vector of length 3, with u/ra for vegetation layers
* if overstory == true the first value is u/ra in vic vegetation,
* and the second value the u/ra in the overstory. otherwise the first
* value is u/ra in the vic vegetation and the second value
* is not used. the third value is u/ra over snow.

      real u(3)           ! wind for vegetation layers (m/s)
      real ra(3)          ! aerodynamic resistance (s/m)
* ----------------------------------------------------------------------

* ------------------- common block variables ---------------------------
#include <vic_veglib.h>
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      integer ost         ! 1 with overstory  0 no overstory
      real z0rl           ! roughness length (m)
      real d              ! vegetation displacement (m)
      real height         ! height of the veg layers (top layer first) (m)
      real trunk          ! multiplier for height that indictaes the tunk

      real refh, d_lower, d_upper, k2
      real uh, ut, uw, z0sl, z0sn
      real z0_lower, z0_upper, zt, zw, tmp_wind, rati

      integer i

      real vnwd           ! attenuation coeff for wind in the overstory
      real von_k          ! von constant
      real huge_resist
      parameter (vnwd=0.5,von_k = 0.41, huge_resist = 1.e2)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* set parameter values
* ----------------------------------------------------------------------

      z0sl  = slrough
      z0sn  = snrough

      k2    = von_k * von_k
      z0rl  = veg_rough(month,vtype)

* z0rl is very important to determine the aerodynamic resistance
* however, in according with vic_veglib.h, the value of veg_rough
* is too big for the some vegetations (6/2004,ji)

      z0rl = z0rl/5.0     ! 6/2004

      ost   = veg_ost(vtype)
      trunk = veg_trunk(vtype)
      d     = veg_d(month,vtype)
      height= d/0.67

* ----------------------------------------------------------------------
* no overstory, thus maximum one soil layer
* ----------------------------------------------------------------------

      if(d.lt.zref) then
         refh = zref
      else
         refh = d + zref + z0rl
         wind = wind*log(refh)/log(zref)
      end if

      do i = 1, 3         ! initialize u and ra
         u(i)  = 0
         ra(i) = 0
      end do

      if (ost.eq.0) then
         z0_lower = z0rl
         d_lower  = d
    
* ----------------------------------------------------------------------
* without snow 
* ----------------------------------------------------------------------
         u(1)  = log((2. + z0_lower)/z0_lower)/
     &           log((refh - d_lower)/z0_lower)

* original vic
c         ra(1) = log((2.+(1.0/0.63-1.0)*d_lower)/z0_lower)
c     &     * log((2.+(1.0/0.63-1.0)*d_lower)/(0.1*z0_lower))/k2
 
         ra(1) = log((2.+(1.0/0.63-1.0)*d_lower)/z0_lower)
     &     * log((refh+(1.0/0.63-1.0)*d_lower)/z0_lower)/k2

* ----------------------------------------------------------------------
* if snow-covered ground (if currently no snow, we still need compute 
* u(3) and ra(3), because there may have snowfalling onto no snow-ground.)
* ----------------------------------------------------------------------
         u(3)  = log((2.+z0sn)/z0sn)/log(refh/z0sn)
         ra(3) = log((2.+z0sn)/z0sn)*log(refh/z0sn)/k2

      else

* ----------------------------------------------------------------------
* overstory present, one or two vegetation layers possible 
* ----------------------------------------------------------------------
         z0_upper = z0rl
         d_upper  = d
    
         z0_lower = z0sl
         d_lower  = 0
    
         zw = 1.5 * height - 0.5 * d_upper
         zt = trunk * height
         if (zt .lt. (z0_lower+d_lower)) then
            print *, 'error: trunk space height below lower boundary'
            call stoprun
         end if

* ----------------------------------------------------------------------
* resistance for overstory
* ----------------------------------------------------------------------

         ra(2) = log((refh-d_upper)/z0_upper)/k2
     &           * (height/(vnwd*(zw-d_upper)) 
     &           * (exp(vnwd*(1-(d_upper+z0_upper)/height))-1)
     &           + (zw-height)/(zw-d_upper)
     &           + log((refh-d_upper)/(zw-d_upper)))
    
* ----------------------------------------------------------------------
* wind at different levels in the profile 
* ----------------------------------------------------------------------

         uw=log((zw-d_upper)/z0_upper)/log((refh-d_upper)/z0_upper)
         uh=uw - (1-(height-d_upper)/(zw-d_upper))
     &      / log((refh-d_upper)/z0_upper)
         u(2) = uh * exp(vnwd * ((z0_upper+d_upper)/height - 1.))
         ut   = uh * exp(vnwd * (zt/height - 1.))
    
* ----------------------------------------------------------------------
* resistance at the lower boundary 
* ----------------------------------------------------------------------

         u(1)  = log((2.+z0_upper)/z0_upper)/log((refh-d_upper)
     &           /z0_upper)

c         ra(1) = log((2. + (1.0/0.63 - 1.0) * d_upper) / z0_upper)
c     &           *log((2.+(1.0/0.63-1.0)*d_upper)/(0.1*z0_upper))/k2

         ra(1) = log((2. + (1.0/0.63 - 1.0) * d_upper) / z0_upper)
     &           *log((refh+(1.0/0.63-1.0)*d_upper)/z0_upper)/k2

* ----------------------------------------------------------------------
* if snow covered ground
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* case 1: the wind profile to a height of 2m above the lower 
*         boundary is entirely logarithmic
* ----------------------------------------------------------------------

         if (zt .gt. (2. + z0sn)) then
            u(3) =ut*log((2.+z0sn)/z0sn)/log(zt/z0sn)
            ra(3)=log((2.+z0sn)/z0sn)*log(zt/z0sn)/(k2*ut) 
               
* ----------------------------------------------------------------------
* case 2: the wind profile to a height of 2m above the lower boundary
*         is part logarithmic and part exponential, but the top of
*         the overstory is more than 2 m above the lower boundary
* ----------------------------------------------------------------------

         else if (height .gt. (2. + z0sn)) then
            u(3)  = uh * exp(vnwd * ((2. + z0sn)/height - 1.))
            ra(3) = log(zt/z0sn) * log(zt/z0sn)/(k2*ut) +
     &           height * log((refh-d_upper)/z0_upper) / 
     &           (vnwd*k2*(zw-d_upper)) *(exp(vnwd*(1-zt/height))
     &           - exp(vnwd*(1-(z0sn+2.)/height)))
            
* ----------------------------------------------------------------------
* case 3: the top of the overstory is less than 2 m above the lower 
*         boundary.  the wind profile above the lower boundary is 
*         part logarithmic and part exponential, but only extends to 
*         the top of the overstory
* ----------------------------------------------------------------------

         else 
            u(3)  = uh
            ra(3) = log(zt/z0sn) * log(zt/z0sn)/(k2*ut) +
     &           height*log((refh-d_upper)/z0_upper) / 
     &           (vnwd*k2*(zw-d_upper))*(exp(vnwd*(1-zt/height))-1)
         end if
      end if


* -- to refrain the aerodynamic resistance, ra=ra*log(ref/2.0) 6/2004 ji

      tmp_wind = max(1.0,min(10.0,wind))
      rati = log(refh/2.0)
      do i = 1, 3
         ra(i) = ra(i)*rati
      enddo

      if(wind.gt.0.) then
         u(1) = u(1)*wind
         ra(1)= ra(1)/tmp_wind

         if(ost.ne.0) then
            u(2) = u(2)*wind
            ra(2)= ra(2)/tmp_wind
         end if
         u(3) = u(3)*wind
         ra(3)= ra(3)/tmp_wind
      else
         u(1) = u(1)*wind
         ra(1) = huge_resist
         if(ost.ne.0) then
            u(2) = u(2)*wind
            ra(2)= huge_resist
         end if
         u(3) = u(3)*wind
         ra(3)= huge_resist
      end if

#ifdef DBGVIC
      print *,'in aerodynamic ra ',ra,' u ',u
      print *,'   aerodynamic ',vtype, month, slrough, snrough,
     &                          zref,    wind
#endif

      return
      end
