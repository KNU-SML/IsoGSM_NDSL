      function snksrc (tavg,smc,sh2o,zsoil,nsoil,
     &                 smcmax,psisat,bexp,dt,k,qtot) 
c      
      implicit none
c      
c ----------------------------------------------------------------------
c function snksrc
c ----------------------------------------------------------------------
c calculate sink/source term of the termal diffusion equation. (sh2o) is
c available liqued water.
c ----------------------------------------------------------------------
      integer k
      integer nsoil
c      
      real bexp
      real df
      real dh2o
      real dt
      real dz
      real dzh
      real free
      real frh2o
      real hlice
      real psisat
      real qtot
      real sh2o
      real smc
      real smcmax
      real snksrc
      real t0
      real tavg
      real tdn
      real tm
      real tup
      real tz
      real x0
      real xdn
      real xh2o
      real xup
      real zsoil (nsoil)
c
      parameter(dh2o = 1.0000e3)
      parameter(hlice = 3.3350e5)
      parameter(t0 = 2.7315e2)
c      
      if (k .eq. 1) then
        dz = -zsoil(1)
      else
        dz = zsoil(k-1)-zsoil(k)
      endif
c ----------------------------------------------------------------------
c via function frh2o, compute potential or 'equilibrium' unfrozen
c supercooled free water for given soil type and soil layer temperature.
c function frh20 invokes eqn (17) from v. koren et al (1999, jgr, vol.
c 104, pg 19573).  (aside:  latter eqn in journal in centigrade units.
c routine frh2o use form of eqn in kelvin units.)
c ----------------------------------------------------------------------
      free = frh2o(tavg,smc,sh2o,smcmax,bexp,psisat)
c ----------------------------------------------------------------------
c in next block of code, invoke eqn 18 of v. koren et al (1999, jgr,
c vol. 104, pg 19573.)  that is, first estimate the new amountof liquid
c water, 'xh2o', implied by the sum of (1) the liquid water at the begin
c of current time step, and (2) the freeze of thaw change in liquid
c water implied by the heat flux 'qtot' passed in from routine hrt.
c second, determine if xh2o needs to be bounded by 'free' (equil amt) or
c if 'free' needs to be bounded by xh2o.
c ----------------------------------------------------------------------
      xh2o = sh2o + qtot*dt/(dh2o*hlice*dz)
c ----------------------------------------------------------------------
c first, if freezing and remaining liquid less than lower bound, then
c reduce extent of freezing, thereby letting some or all of heat flux
c qtot cool the soil temp later in routine hrt.
c ----------------------------------------------------------------------
      if ( xh2o .lt. sh2o .and. xh2o .lt. free) then 
        if ( free .gt. sh2o ) then
          xh2o = sh2o
        else
          xh2o = free
        endif
      endif
c ----------------------------------------------------------------------
c second, if thawing and the increase in liquid water greater than upper
c bound, then reduce extent of thaw, thereby letting some or all of heat
c flux qtot warm the soil temp later in routine hrt.
c ----------------------------------------------------------------------
      if ( xh2o .gt. sh2o .and. xh2o .gt. free )  then
        if ( free .lt. sh2o ) then
          xh2o = sh2o
        else
          xh2o = free
        endif
      endif 
c
      if (xh2o .lt. 0.) xh2o = 0.
      if (xh2o .gt. smc) xh2o = smc
c ----------------------------------------------------------------------
c calculate phase-change heat source/sink term for use in routine hrt
c and update liquid water to reflcet final freeze/thaw increment.
c ----------------------------------------------------------------------
      snksrc = -dh2o*hlice*dz*(xh2o-sh2o)/dt
      sh2o = xh2o
c ----------------------------------------------------------------------
c end function snksrc
c ----------------------------------------------------------------------
77    return
      end
