      subroutine canres (solar,ch,sfctmp,q2,sfcprs,smc,zsoil,nsoil,
     &                   smcwlt,smcref,rsmin,rc,pc,nroot,q2sat,dqsdt2, 
     &                   topt,rsmax,rgl,hs,xlai,
     &                   rcs,rct,rcq,rcsoil)
c
      implicit none
c ----------------------------------------------------------------------
c subroutine canres                    
c ----------------------------------------------------------------------
c calculate canopy resistance which depends on incoming solar radiation,
c air temperature, atmospheric water vapor pressure deficit at the
c lowest model level, and soil moisture (preferably unfrozen soil
c moisture rather than total)
c ----------------------------------------------------------------------
c source:  jarvis (1976), noilhan and planton (1989, mwr), jacquemin and
c noilhan (1990, blm)
c see also:  chen et al (1996, jgr, vol 101(d3), 7251-7268), eqns 12-14
c and table 2 of sec. 3.1.2         
c ----------------------------------------------------------------------
c input:
c   solar   incoming solar radiation
c   ch      surface exchange coefficient for heat and moisture
c   sfctmp  air temperature at 1st level above ground
c   q2      air humidity at 1st level above ground
c   q2sat   saturation air humidity at 1st level above ground
c   dqsdt2  slope of saturation humidity function wrt temp
c   sfcprs  surface pressure
c   smc     volumetric soil moisture 
c   zsoil   soil depth (negative sign, as it is below ground)
c   nsoil   no. of soil layers
c   nroot   no. of soil layers in root zone (1.le.nroot.le.nsoil)
c   xlai    leaf area index
c   smcwlt  wilting point
c   smcref  reference soil moisture (where soil water deficit stress
c             sets in)
c rsmin, rsmax, topt, rgl, hs are canopy stress parameters set in
c   surboutine redprm
c output:
c   pc  plant coefficient
c   rc  canopy resistance
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer k
      integer nroot
      integer nsoil
c
      real ch
      real cp
      real delta
      real dqsdt2
      real ff
      real gx
      real hs
      real p
      real part(nsold) 
      real pc
      real q2
      real q2sat
      real rc
      real rsmin
      real rcq
      real rcs
      real rcsoil
      real rct
      real rd
      real rgl
      real rr
      real rsmax
      real sfcprs
      real sfctmp
      real sigma
      real slv
      real smc(nsoil)
      real smcref
      real smcwlt
      real solar
      real topt
      real slvcp
      real st1
      real tair4
      real xlai
      real zsoil(nsoil)
c
      parameter(cp = 1004.5)
      parameter(rd = 287.04)
      parameter(sigma = 5.67e-8)
      parameter(slv=2.501000e6)
c ----------------------------------------------------------------------
c initialize canopy resistance multiplier terms.
c ----------------------------------------------------------------------
      rcs = 0.0
      rct = 0.0
      rcq = 0.0
      rcsoil = 0.0
      rc = 0.0
c ----------------------------------------------------------------------
c contribution due to incoming solar radiation
c ----------------------------------------------------------------------
      ff = 0.55*2.0*solar/(rgl*xlai)
      rcs = (ff + rsmin/rsmax) / (1.0 + ff)
      rcs = max(rcs,0.0001)
c ----------------------------------------------------------------------
c contribution due to air temperature at first model level above ground
c rct expression from noilhan and planton (1989, mwr).
c ----------------------------------------------------------------------
      rct = 1.0 - 0.0016*((topt-sfctmp)**2.0)
      rct = max(rct,0.0001)
c ----------------------------------------------------------------------
c contribution due to vapor pressure deficit at first model level.
c rcq expression from ssib 
c ----------------------------------------------------------------------
      rcq = 1.0/(1.0+hs*(q2sat-q2))
      rcq = max(rcq,0.01)
c ----------------------------------------------------------------------
c contribution due to soil moisture availability.
c determine contribution from each soil layer, then add them up.
c ----------------------------------------------------------------------
      gx = (smc(1) - smcwlt) / (smcref - smcwlt)
      if (gx .gt. 1.) gx = 1.
      if (gx .lt. 0.) gx = 0.
c ----------------------------------------------------------------------
c use soil depth as weighting factor
c ----------------------------------------------------------------------
      part(1) = (zsoil(1)/zsoil(nroot)) * gx
c ----------------------------------------------------------------------
c use root distribution as weighting factor
c      part(1) = rtdis(1) * gx
c ----------------------------------------------------------------------
      do k = 2,nroot
        gx = (smc(k) - smcwlt) / (smcref - smcwlt)
        if (gx .gt. 1.) gx = 1.
        if (gx .lt. 0.) gx = 0.
c ----------------------------------------------------------------------
c use soil depth as weighting factor        
c ----------------------------------------------------------------------
        part(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot)) * gx
c ----------------------------------------------------------------------
c use root distribution as weighting factor
c        part(k) = rtdis(k) * gx 
c ----------------------------------------------------------------------
      end do
      do k = 1,nroot
        rcsoil = rcsoil+part(k)
      end do
      rcsoil = max(rcsoil,0.0001)
c ----------------------------------------------------------------------
c determine canopy resistance due to all factors.  convert canopy
c resistance (rc) to plant coefficient (pc) to be used with potential
c evap in determining actual evap.  pc is determined by:
c   pc * linerized penman potential evap =
c   penman-monteith actual evaporation (containing rc term).
c ----------------------------------------------------------------------
      rc = rsmin/(xlai*rcs*rct*rcq*rcsoil)
c      tair4 = sfctmp**4.
c      st1 = (4.*sigma*rd)/cp
c      slvcp = slv/cp
c      rr = st1*tair4/(sfcprs*ch) + 1.0
      rr = (4.*sigma*rd/cp)*(sfctmp**4.)/(sfcprs*ch) + 1.0
c-mk
      delta = (slv/cp)*dqsdt2
      pc = (rr+delta)/(rr*(1.+rc*ch)+delta)
c ----------------------------------------------------------------------
c end subroutine canres
c ----------------------------------------------------------------------
      return
      end
