      subroutine sstep (sh2oout,sh2oin,cmc,rhstt,rhsct,dt,
     &                  nsoil,smcmax,cmcmax,runoff3,zsoil,smc,sice,
     &                  ai,bi,ci)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine sstep
c ----------------------------------------------------------------------
c calculate/update soil moisture content values and canopy moisture
c content values.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer i
      integer k 
      integer kk11
      integer nsoil
c
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
      real ciin(nsold)
      real cmc
      real cmcmax
      real ddz
      real dt
      real rhsct
      real rhstt(nsoil)
      real rhsttin(nsoil)
      real runoff3
      real sh2oin(nsoil)
      real sh2oout(nsoil)
      real sice(nsoil)
      real smc(nsoil)
      real smcmax
      real stot
      real wplus
      real zsoil(nsoil)
c ----------------------------------------------------------------------
c create 'amount' values of variables to be input to the
c tri-diagonal matrix routine.
c ----------------------------------------------------------------------
      do k = 1,nsoil
        rhstt(k) = rhstt(k) * dt
        ai(k) = ai(k) * dt
        bi(k) = 1. + bi(k) * dt
        ci(k) = ci(k) * dt
      end do
c ----------------------------------------------------------------------
c copy values for input variables before call to rosr12
c ----------------------------------------------------------------------
      do k = 1,nsoil
        rhsttin(k) = rhstt(k)
      end do
      do k = 1,nsold
        ciin(k) = ci(k)
      end do
c ----------------------------------------------------------------------
c call rosr12 to solve the tri-diagonal matrix
c ----------------------------------------------------------------------
      call rosr12 (ci,ai,bi,ciin,rhsttin,rhstt,nsoil)
c ----------------------------------------------------------------------
c sum the previous smc value and the matrix solution to get a
c new value.  min allowable value of smc will be 0.02.
c runoff3: runoff within soil layers
c ----------------------------------------------------------------------
      wplus = 0.0
      runoff3 = 0.
      ddz = -zsoil(1)
      do k = 1,nsoil
        if (k .ne. 1) ddz = zsoil(k - 1) - zsoil(k)
        sh2oout(k) = sh2oin(k) + ci(k) + wplus / ddz
c
        stot = sh2oout(k) + sice(k)
        if (stot .gt. smcmax) then
          if (k .eq. 1) then
            ddz = -zsoil(1)
          else
            kk11 = k - 1
            ddz = -zsoil(k) + zsoil(kk11)
          endif
          wplus = (stot-smcmax) * ddz
        else
          wplus = 0.
        endif
        smc(k) = max ( min(stot,smcmax),0.02 )
        sh2oout(k) = max((smc(k)-sice(k)),0.0)
      end do
      runoff3 = wplus
c ----------------------------------------------------------------------
c update canopy water content/interception (cmc).  convert rhsct to 
c an 'amount' value and add to previous cmc value to get new cmc.
c ----------------------------------------------------------------------
      cmc = cmc + dt * rhsct
      if (cmc .lt. 1.e-20) cmc=0.0
      cmc = min(cmc,cmcmax)
c ----------------------------------------------------------------------
c end subroutine sstep
c ----------------------------------------------------------------------
      return
      end
