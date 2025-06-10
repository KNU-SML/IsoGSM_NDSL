      subroutine transp (et,nsoil,etp1,smc,cmc,zsoil,shdfac,smcwlt,
     &                   cmcmax,pc,cfactr,smcref,sfctmp,q2,nroot,rtdis)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine transp
c ----------------------------------------------------------------------
c calculate transpiration for the veg class.
c ----------------------------------------------------------------------
      integer i
      integer k
      integer nsoil
      integer nroot
c
      real cfactr
      real cmc
      real cmcmax
      real denom
      real et(nsoil)
      real etp1
      real etp1a
      real gx (7)
c.....real part(nsoil)
      real pc
      real q2
      real rtdis(nsoil)
      real rtx
      real sfctmp
      real sgx
      real shdfac
      real smc(nsoil)
      real smcref
      real smcwlt
      real zsoil(nsoil)
c ----------------------------------------------------------------------
c initialize plant transp to zero for all soil layers.
c ----------------------------------------------------------------------
      do k = 1,nsoil
        et(k) = 0.
      end do
c ----------------------------------------------------------------------
c calculate an 'adjusted' potential transpiration
c if statement below to avoid tangent linear problems near zero
c note: gx and other terms below redistribute transpiration by layer,
c et(k), as a function of soil moisture availability, while preserving
c total etp1a.
c ----------------------------------------------------------------------
      if (cmc .ne. 0.0) then
        etp1a = shdfac * pc * etp1 * (1.0 - (cmc /cmcmax) ** cfactr)
      else
        etp1a = shdfac * pc * etp1
      endif
      sgx = 0.0
      do i = 1,nroot
        gx(i) = ( smc(i) - smcwlt ) / ( smcref - smcwlt )
        gx(i) = max ( min ( gx(i), 1. ), 0. )
        sgx = sgx + gx (i)
      end do
      sgx = sgx / nroot
      denom = 0.
      do i = 1,nroot
        rtx = rtdis(i) + gx(i) - sgx
        gx(i) = gx(i) * max ( rtx, 0. )
        denom = denom + gx(i)
      end do
      if (denom .le. 0.0) denom = 1.
      do i = 1,nroot
        et(i) = etp1a * gx(i) / denom
      end do
c ----------------------------------------------------------------------
c above code assumes a vertically uniform root distribution
c code below tests a variable root distribution
c ----------------------------------------------------------------------
c      et(1) = ( zsoil(1) / zsoil(nroot) ) * gx * etp1a
c      et(1) = ( zsoil(1) / zsoil(nroot) ) * etp1a
c ----------------------------------------------------------------------
c using root distribution as weighting factor
c ----------------------------------------------------------------------
c      et(1) = rtdis(1) * etp1a
c      et(1) = etp1a * part(1)
c ----------------------------------------------------------------------
c loop down thru the soil layers repeating the operation above,
c but using the thickness of the soil layer (rather than the
c absolute depth of each layer) in the final calculation.
c ----------------------------------------------------------------------
c      do k = 2,nroot
c        gx = ( smc(k) - smcwlt ) / ( smcref - smcwlt )
c        gx = max ( min ( gx, 1. ), 0. )
c test canopy resistance
c        gx = 1.0
c        et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*gx*etp1a
c        et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*etp1a
c ----------------------------------------------------------------------
c using root distribution as weighting factor
c ----------------------------------------------------------------------
c        et(k) = rtdis(k) * etp1a
c        et(k) = etp1a*part(k)
c      end do      
c ----------------------------------------------------------------------
c end subroutine transp
c ----------------------------------------------------------------------
      return
      end
