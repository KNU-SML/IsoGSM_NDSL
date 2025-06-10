      function frh2o (tkelv,smc,sh2o,smcmax,bexp,psis)
      implicit none
c ----------------------------------------------------------------------
c function frh2o
c ----------------------------------------------------------------------
c calculate amount of supercooled liquid soil water content if
c temperature is below 273.15k (t0).  requires newton-type iteration to
c solve the nonlinear implicit equation given in eqn 17 of koren et al
c (1999, jgr, vol 104(d16), 19569-19585).
c ----------------------------------------------------------------------
c new version (june 2001): much faster and more accurate newton
c iteration achieved by first taking log of eqn cited above -- less than
c 4 (typically 1 or 2) iterations achieves convergence.  also, explicit
c 1-step solution option for special case of parameter ck=0, which
c reduces the original implicit equation to a simpler explicit form,
c known as the "flerchinger eqn". improved handling of solution in the
c limit of freezing point temperature t0.
c ----------------------------------------------------------------------
c input:
c
c   tkelv.........temperature (kelvin)
c   smc...........total soil moisture content (volumetric)
c   sh2o..........liquid soil moisture content (volumetric)
c   smcmax........saturation soil moisture content (from redprm)
c   b.............soil type "b" parameter (from redprm)
c   psis..........saturated soil matric potential (from redprm)
c
c output:
c   frh2o.........supercooled liquid water content
c ----------------------------------------------------------------------
      real bexp
      real blim
      real bx
      real ck
      real denom
      real df
      real dh2o
      real dice
      real dswl
      real error
      real fk
      real frh2o
      real gs
      real hlice
      real psis
      real sh2o
      real smc
      real smcmax
      real swl
      real swlk
      real tkelv
      real t0
      integer nlog
      integer kcount
      parameter(ck = 8.0)
c      parameter(ck = 0.0)
      parameter(blim = 5.5)
      parameter(error = 0.005)
      parameter(hlice = 3.335e5)
      parameter(gs = 9.81)
      parameter(dice = 920.0)
      parameter(dh2o = 1000.0)
      parameter(t0 = 273.15)
c ----------------------------------------------------------------------
c limits on parameter b: b < 5.5  (use parameter blim)
c simulations showed if b > 5.5 unfrozen water content is
c non-realistically high at very low temperatures.
c ----------------------------------------------------------------------
      bx = bexp
      if (bexp .gt. blim) bx = blim
c ----------------------------------------------------------------------
c initializing iterations counter and iterative solution flag.
c ----------------------------------------------------------------------
      nlog=0
      kcount=0
c ----------------------------------------------------------------------
c  if temperature not significantly below freezing (t0), sh2o = smc
c ----------------------------------------------------------------------
      if (tkelv .gt. (t0 - 1.e-3)) then
        frh2o = smc
      else
        if (ck .ne. 0.0) then
c ----------------------------------------------------------------------
c option 1: iterated solution for nonzero ck
c in koren et al, jgr, 1999, eqn 17
c ----------------------------------------------------------------------
c initial guess for swl (frozen content)
c ----------------------------------------------------------------------
          swl = smc-sh2o
c ----------------------------------------------------------------------
c keep within bounds.
c ----------------------------------------------------------------------
          if (swl .gt. (smc-0.02)) swl = smc-0.02
          if (swl .lt. 0.) swl = 0.
c ----------------------------------------------------------------------
c  start of iterations
c ----------------------------------------------------------------------
          do while ( (nlog .lt. 10) .and. (kcount .eq. 0) )
            nlog = nlog+1
            df = alog(( psis*gs/hlice ) * ( ( 1.+ck*swl )**2. ) *
     &        ( smcmax/(smc-swl) )**bx) - alog(-(tkelv-t0)/tkelv)
            denom = 2. * ck / ( 1.+ck*swl ) + bx / ( smc - swl )
            swlk = swl - df/denom
c ----------------------------------------------------------------------
c bounds useful for mathematical solution.
c ----------------------------------------------------------------------
            if (swlk .gt. (smc-0.02)) swlk = smc - 0.02
            if (swlk .lt. 0.) swlk = 0.
c ----------------------------------------------------------------------
c mathematical solution bounds applied.
c ----------------------------------------------------------------------
            dswl = abs(swlk-swl)
            swl = swlk
c ----------------------------------------------------------------------
c if more than 10 iterations, use explicit method (ck=0 approx.)
c when dswl less or eq. error, no more iterations required.
c ----------------------------------------------------------------------
            if ( dswl .le. error )  then
              kcount = kcount+1
            endif
          end do
c ----------------------------------------------------------------------
c  end of iterations
c ----------------------------------------------------------------------
c bounds applied within do-block are valid for physical solution.
c ----------------------------------------------------------------------
          frh2o = smc - swl
c ----------------------------------------------------------------------
c end option 1
c ----------------------------------------------------------------------
        endif
c ----------------------------------------------------------------------
c option 2: explicit solution for flerchinger eq. i.e. ck=0
c in koren et al., jgr, 1999, eqn 17
c apply physical bounds to flerchinger solution
c ----------------------------------------------------------------------
        if (kcount .eq. 0) then
clu.....comment out the following line to shorten the standard output
!clu      print*,'flerchinger used in new version. iterations=',nlog
          fk = (((hlice/(gs*(-psis)))*
     &      ((tkelv-t0)/tkelv))**(-1/bx))*smcmax
          if (fk .lt. 0.02) fk = 0.02
          frh2o = min (fk, smc)
c ----------------------------------------------------------------------
c end option 2
c ----------------------------------------------------------------------
        endif
      endif
c ----------------------------------------------------------------------
c end function frh2o
c ----------------------------------------------------------------------
      return
      end
