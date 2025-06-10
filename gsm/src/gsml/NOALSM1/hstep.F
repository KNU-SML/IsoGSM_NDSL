      subroutine hstep (stcout,stcin,rhsts,dt,nsoil,ai,bi,ci)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine hstep
c ----------------------------------------------------------------------
c calculate/update the soil temperature field.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer k
      integer nsoil
c
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
      real ciin(nsold)
      real dt
      real rhsts(nsoil)
      real rhstsin(nsoil)
      real stcin(nsoil)
      real stcout(nsoil)
c ----------------------------------------------------------------------
c create finite difference values for use in rosr12 routine
c ----------------------------------------------------------------------
      do k = 1,nsoil
        rhsts(k) = rhsts(k) * dt
        ai(k) = ai(k) * dt
        bi(k) = 1. + bi(k) * dt
        ci(k) = ci(k) * dt
      end do
c ----------------------------------------------------------------------
c copy values for input variables before call to rosr12
c ----------------------------------------------------------------------
      do k = 1,nsoil
         rhstsin(k) = rhsts(k)
      end do
      do k = 1,nsold
        ciin(k) = ci(k)
      end do
c ----------------------------------------------------------------------
c solve the tri-diagonal matrix equation
c ----------------------------------------------------------------------
      call rosr12(ci,ai,bi,ciin,rhstsin,rhsts,nsoil)
c ----------------------------------------------------------------------
c calc/update the soil temps using matrix solution
c ----------------------------------------------------------------------
      do k = 1,nsoil
        stcout(k) = stcin(k) + ci(k)
      end do
c ----------------------------------------------------------------------
c end subroutine hstep
c ----------------------------------------------------------------------
      return
      end
