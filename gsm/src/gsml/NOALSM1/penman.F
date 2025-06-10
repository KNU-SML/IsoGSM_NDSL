      subroutine penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,t24,ssoil,
     &                   q2,q2sat,etp,rch,epsca,rr,snowng,frzgra,
     &                   dqsdt2,flx2)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine penman
c ----------------------------------------------------------------------
c calculate potential evaporation for the current point.  various
c partial sums/products are also calculated and passed back to the
c calling routine for later use.
c ----------------------------------------------------------------------
      logical snowng
      logical frzgra
c
      real a
      real beta
      real ch
      real cp
      real cph2o
      real cpice
      real delta
      real dqsdt2
      real elcp
      real epsca
      real etp
      real fdown
      real flx2
      real fnet
      real lsubc
      real lsubf
      real prcp
      real q2
      real q2sat
      real r
      real rad
      real rch
      real rho
      real rr
      real ssoil
      real sfcprs
      real sfctmp
      real sigma
      real t24
      real t2v
      real th2
c
      parameter(cp = 1004.6)
      parameter(cph2o = 4.218e+3)
      parameter(cpice = 2.106e+3)
      parameter(r = 287.04)
      parameter(elcp = 2.4888e+3)
      parameter(lsubf = 3.335e+5)
      parameter(lsubc = 2.501000e+6)
      parameter(sigma = 5.67e-8)
c ----------------------------------------------------------------------
c executable code begins here:
c ----------------------------------------------------------------------
      flx2 = 0.0
c ----------------------------------------------------------------------
c prepare partial quantities for penman equation.
c ----------------------------------------------------------------------
      delta = elcp * dqsdt2
      t24 = sfctmp * sfctmp * sfctmp * sfctmp
      rr = t24 * 6.48e-8 /(sfcprs * ch) + 1.0
      rho = sfcprs / (r * t2v)
      rch = rho * cp * ch
c ----------------------------------------------------------------------
c adjust the partial sums / products with the latent heat
c effects caused by falling precipitation.
c ----------------------------------------------------------------------
      if (.not. snowng) then
        if (prcp .gt. 0.0) rr = rr + cph2o*prcp/rch
      else
        rr = rr + cpice*prcp/rch
      endif
      fnet = fdown - sigma*t24 - ssoil
c ----------------------------------------------------------------------
c include the latent heat effects of frzng rain converting to ice on
c impact in the calculation of flx2 and fnet.
c ----------------------------------------------------------------------
      if (frzgra) then
        flx2 = -lsubf * prcp
        fnet = fnet - flx2
      endif
c ----------------------------------------------------------------------
c finish penman equation calculations.
c ----------------------------------------------------------------------
      rad = fnet/rch + th2 - sfctmp
      a = elcp * (q2sat - q2)
      epsca = (a*rr + rad*delta) / (delta + rr)
      etp = epsca * rch / lsubc
c ----------------------------------------------------------------------
c end subroutine penman
c ----------------------------------------------------------------------
      return
      end
