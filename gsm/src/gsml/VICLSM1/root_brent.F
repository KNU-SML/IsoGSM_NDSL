      function root_brent(ubound, lbound, ctrfct)

#include <vartyp.h>

************************************************************************
* source file:       root_brent.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (root_brent.c)
*
* abstract: this function returns the temperature for 
*  which the sum of the energy balance terms is zero in the 
*  interval [mint, maxt]. the temperature is calculated 
*  to within a tolerance (6 * macheps * |t| + 2 * t), where macheps 
*  is the relative machine precision and t is a positive tolerance, as 
*  specified.
*            the function assures that f(mintsurf) and f(maxtsurf) have 
*  opposite signs. if this is not the case the program will stop.  in 
*  addition the program will perform not more than a certain number of 
*  iterations, as specified in brent.h, and will abort if more iterations 
*  are needed.
*                 
*  general documentation for this module
*
*  source: brent, r. p., 1973, algorithms for minimization without
*   derivatives,prentice hall, inc., englewood cliffs, new jersey chapter 4
*
*  the method is also discussed in:
*  press, w. h., s. a. teukolsky, w. t. vetterling, b. p. flannery, 1992,
*  numerical recipes in fortran, the art of scientific computing,
*  second edition, cambridge university press
*  (be aware that this book discusses a brent method for minimization (brent), 
*  and one for root finding (zbrent).  the latter one is similar to the one 
*  implemented here and is also copied from brent [1973].)
************************************************************************

* ------------------- input variables ----------------------------------
      real ubound           ! upper bound for root
      real lbound           ! lower bound for root

      character*(*)  ctrfct ! determine function
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      real a, b, c, d, e, fa, fb, fc
      real m, p, q, r, s, tol

      real root_brent         ! function name
      real func_surfeb        ! function name
      real func_soiltemp      ! function name
      real func_snoweb        ! function name

      real tstep, tt
      parameter (tstep = 0.5, tt = 1.e-7)

      integer maxtries, maxiter
      parameter (maxtries=100, maxiter=1000)
      real macheps
      parameter (macheps=3.e-8)

      integer i, j            ! loop index

* ----------------------------------------------------------------------
* initialize variable argujment list
* ----------------------------------------------------------------------

      a = lbound
      b = ubound
      
      if(ctrfct.eq.'func_surfeb')then
         fa = func_surfeb(a)
         fb = func_surfeb(b)
      elseif(ctrfct.eq.'func_soiltemp')then
         fa = func_soiltemp(a)
         fb = func_soiltemp(b)
      elseif(ctrfct.eq.'func_snoweb')then
         fa = func_snoweb(a)
         fb = func_snoweb(b)
      end if

* ----------------------------------------------------------------------
* if root not bracketed attempt to bracket the root
* ----------------------------------------------------------------------

      j = 0

      do while ((fa * fb).ge. 0  .and. j .lt. maxtries)

         a = a-tstep
         if(ctrfct.eq.'func_snoweb') then
            b = b
         else
            b = b+tstep
         endif

         if(ctrfct.eq.'func_surfeb')then
            fa = func_surfeb(a)
            fb = func_surfeb(b)
         elseif(ctrfct.eq.'func_soiltemp')then
            fa = func_soiltemp(a)
            fb = func_soiltemp(b)
         elseif(ctrfct.eq.'func_snoweb')then
            fa = func_snoweb(a)
         end if

         j = j + 1
      end do

      if ((fa * fb) .ge. 0) then
         print *, 'error: first error in root_brent fa * fb >= 0'
         print *, 'ctrfct ',ctrfct, ' a=',a,' b=',b
         print *, ' fa=',fa,' fb=',fb
         call stoprun
      end if
  
      fc = fb
      
      do i = 1, maxiter

         if(fb*fc .gt. 0) then
            c = a
            fc = fa
            d = b - a
            e = d
         end if
    
         if (abs(fc) .lt. abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         end if
      
         tol = 2 * macheps * abs(b) + tt
         m = 0.5 * (c - b)
    
         if (abs(m) .le. tol .or. fb.eq.0.0) then
            if(ctrfct.eq.'func_surfeb')then
               fb = func_surfeb(b)
            elseif(ctrfct.eq.'func_soiltemp')then
               fb = func_soiltemp(b)
            elseif(ctrfct.eq.'func_snoweb')then
               fb = func_snoweb(b)
            end if
            root_brent = b
            return
         else
            if (abs(e) .lt. tol .or. abs(fa) .le. abs(fb)) then
               d = m
               e = d
            else 
               s = fb/fa
               if (a.eq.c) then

* ----------------------------------------------------------------------
* linear interpolation           
* ----------------------------------------------------------------------

                  p = 2.0 * m * s
                  q = 1.0 - s
               else

* ----------------------------------------------------------------------
* inverse quadratic interpolation
* ----------------------------------------------------------------------

                  q = fa/fc
                  r = fb/fc
                  p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                  q = (q - 1.0) * (r - 1.0) * (s - 1.0)
               end if
	
               if (p .gt. 0.0) then
                  q = -q
               else
                  p = -p
               end if
               s = e
               e = d
               if(((2.*p).lt.(3.0*m*q-abs(tol*q))).and.
     &            (p.lt.abs(0.5*s*q)))then
                  d = p/q
               else 
                  d = m
                  e = d
               end if
            end if
            a = b
            fa = fb
            if(abs(d).gt.tol) then
               b = b + d
            else
               if(m.gt.0) then
                  b = b + tol
               else
                  b = b - tol
               end if
            end if

            if(ctrfct.eq.'func_surfeb')then
               fb = func_surfeb(b)
            elseif(ctrfct.eq.'func_soiltemp')then
               fb = func_soiltemp(b)
            elseif(ctrfct.eq.'func_snoweb')then
               fb = func_snoweb(b)
            end if

         end if
      end do

      print *, 'error: second error in root_brent too many iterations'
      print *, 'function ctrfct ',ctrfct(1:14), ' a=',a,' b=',b
      print *, 'fa=',fa,' fb=',fb
      print *, 'dumping input variables - check for valid values'
      call stoprun

      return
      end
