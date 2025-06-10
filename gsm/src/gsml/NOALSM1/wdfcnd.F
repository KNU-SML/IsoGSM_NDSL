      subroutine wdfcnd (wdf,wcnd,smc,smcmax,bexp,dksat,dwsat,
     &                   sicemax)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine wdfcnd
c ----------------------------------------------------------------------
c calculate soil water diffusivity and soil hydraulic conductivity.
c ----------------------------------------------------------------------
      real bexp
      real dksat
      real dwsat
      real expon
      real factr1
      real factr2
      real sicemax
      real smc
      real smcmax
      real vkwgt
      real wcnd
      real wdf
c ----------------------------------------------------------------------
c     calc the ratio of the actual to the max psbl soil h2o content
c ----------------------------------------------------------------------
      smc = smc
      smcmax = smcmax
      factr1 = 0.2 / smcmax
      factr2 = smc / smcmax
c ----------------------------------------------------------------------
c prep an expntl coef and calc the soil water diffusivity
c ----------------------------------------------------------------------
      expon = bexp + 2.0
      wdf = dwsat * factr2 ** expon
c ----------------------------------------------------------------------
c frozen soil hydraulic diffusivity.  very sensitive to the vertical
c gradient of unfrozen water. the latter gradient can become very
c extreme in freezing/thawing situations, and given the relatively 
c few and thick soil layers, this gradient sufferes serious 
c trunction errors yielding erroneously high vertical transports of
c unfrozen water in both directions from huge hydraulic diffusivity.  
c therefore, we found we had to arbitrarily constrain wdf 
c --
c version d_10cm: ........  factr1 = 0.2/smcmax
c weighted approach...................... pablo grunmann, 28_sep_1999.
c ----------------------------------------------------------------------
      if (sicemax .gt. 0.0)  then
        vkwgt = 1./(1.+(500.*sicemax)**3.)
        wdf = vkwgt*wdf + (1.- vkwgt)*dwsat*factr1**expon
      endif
c ----------------------------------------------------------------------
c reset the expntl coef and calc the hydraulic conductivity
c ----------------------------------------------------------------------
      expon = (2.0 * bexp) + 3.0
      wcnd = dksat * factr2 ** expon
c ----------------------------------------------------------------------
c end subroutine wdfcnd
c ----------------------------------------------------------------------
      return
      end
