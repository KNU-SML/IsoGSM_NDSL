      subroutine soil_temprfl(nsl,  dtime,  alpnd,  betnd,
     &                      gamnd,     t0,  kapnd,   csnd,
     &                       smnd,  icend,  smxnd,  bubnd,
     &                      expnd,     tn)

#include <vartyp.h>

************************************************************************
* source file:       soil_temprfl.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august, and september  2003
* prgm history:      modified from 4.0.3 uw vic (frozen_soil.c)
*
* abstract: this subroutine iteratively solve the soil temperature
*           profile using a numerical difference equation. the solution
*           equation is second order in space, and first order in time.
*
* comments:
*    please check eqns (3), (4), (5), (6), (7) and (8) in cherkauer &
*    lettenmaier (jgr 1999) for the computing soiltemp profile
************************************************************************

* ------------------- input variables ----------------------------------
      integer nsl           ! number of soil thermal nodes

      real dtime            ! time step (s)
      real alpnd(nsl)       ! coeff for computing soiltemp profile
      real betnd(nsl)       ! coeff for computing soiltemp profile
      real gamnd(nsl)       ! coeff for computing soiltemp profile
      real t0(nsl)          ! olde soil temperature profile (k)
      real kapnd(nsl)       ! soil layer thermal conductivity (w/m/k)
      real csnd(nsl)        ! soil layer heat capacity (j/m^3/k)
      real smnd(nsl)        ! soil smndure at node (mm/mm)
      real icend(nsl)       ! soil icend content at node (mm/mm)
      real smxnd(nsl)       ! maximum soil moisture at node (mm/mm)
      real bubnd(nsl)       ! bubbling pressure of soil (cm)
      real expnd(nsl)       ! parameter for ksat with soil moisture (n/a)   
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real tn(nsl)          ! new soil temperature profile (k)
* ----------------------------------------------------------------------

* ------------------- common block -------------------------------------
#include <vic_soileb.h>
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      integer j, itcount

      real maxdiff, diff, oldt, alpnd2(nsl)
      real a(nsl), b(nsl), c(nsl), d(nsl), e(nsl)

      real root_brent       ! function name

      character*80  ctrfct  ! determine function

      integer maxit
      real soildt           ! used to bracket soil temperatures while
      real threshold        ! temperature profile iteration threshold
      parameter (maxit = 10, soildt=0.25, threshold = 1.e-2)


      real rhoice           ! density of ice (kg/m^3)
      real lf               ! latent heat of freezing (j/kg) at 0c
      real t0c              ! ice/water mix temperature (k)
      parameter (rhoice = 917.0,lf=3.337e5,t0c=273.15)
* ----------------------------------------------------------------------

      do j = 1, nsl
         alpnd2(j) = alpnd(j)*alpnd(j)
      end do
      
      do j = 2, nsl
         if(j.eq.nsl) then
            a(j) = betnd(j)*dtime*(kapnd(j)-kapnd(j-1))
         else
            a(j) = betnd(j)*dtime*(kapnd(j+1)-kapnd(j-1))
         endif
         b(j) = 2.*alpnd2(j)*dtime*kapnd(j)
         c(j) = alpnd2(j)*betnd(j)*csnd(j)*(t0(j)-t0c)
         d(j) = alpnd2(j)*betnd(j)*rhoice*lf
         e(j) = alpnd2(j)*betnd(j)*csnd(j)
     &             + 4.*kapnd(j)*alpnd2(j)*dtime
      end do

* tn(1) = ts, initialized in func_surfeb

      do j = 2, nsl
         tn(j)=t0(j)
      end do

      itcount = 0

 5000 if(itcount.lt.maxit) then
         itcount = itcount + 1
         maxdiff = threshold
         do j= 2, nsl
            oldt = tn(j)
* ----------------------------------------------------------------------
* 2nd order variable kapnd equation 
* ----------------------------------------------------------------------

            if(j.eq.nsl) then
               tl = tn(j)
            else
               tl = tn(j+1)
            endif
            tu = tn(j-1)

            fprime = (tl-tu)/alpnd(j)
            
            if(tn(j).gt.t0c) then
               tn(j) = (a(j)*(tl-tu)
     &                + b(j)*(tl+tu-2*t0c-gamnd(j)*fprime)
     &                + c(j) + d(j)*(0.-icend(j))) / e(j) + t0c
            else
               soilm  = smnd(j)
               soilmx = smxnd(j)
               bubble = bubnd(j)
               expt   = expnd(j)
               ice0   = icend(j)
               gamma  = gamnd(j)
               aa = a(j)
               bb = b(j)
               cc = c(j)
               dd = d(j)
               ee = e(j)
               ctrfct = 'func_soiltemp'

#ifdef DBGVIC
c               print *,'in soil_temprfl bef root_brent',j,aa,bb,cc,dd,ee
#endif

               tn(j) =root_brent(t0(j)+soildt, t0(j)-soildt, ctrfct)
               if(abs(t0(j)-tn(j)).gt.5) then
                  print *,'in soil_temprfl aft ',j,t0(j),tn(j)
c                  stop
               endif
            end if
	
            diff= abs(oldt-tn(j))
            if(diff .gt. maxdiff) maxdiff=diff
         end do
      end if

      if(maxdiff .gt. threshold) then
         goto 5000
      end if

      return
      end
