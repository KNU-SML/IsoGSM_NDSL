      subroutine qsatcn(tt,p,q,dqdt)

!     implicit none
      real tt, p, q, dqdt

      real        airmw, h2omw, lice, runiv, rgas
      parameter ( airmw  = 28.97      )
      parameter ( h2omw  = 18.01      )
      parameter ( lice   = 2.834e6    )
      parameter ( runiv  = 8314.3     )
      parameter ( rgas   = runiv/airmw)
!

      real        esfac, erfac
      parameter ( esfac = h2omw/airmw       )
      parameter ( erfac = (1.0-esfac)/esfac )

      parameter(nx=7501)
      dimension tbsvp(nx), tbdsvp(nx)
      common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp
!
cfpp$ noconcur r
      real d, qx, dqx
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xj  = min(max(c1xsvp+c2xsvp*tt,1.),float(nx))
      jx  = min(xj,nx-1.)
      qx  = tbsvp(jx)  + (xj-jx)*(tbsvp(jx+1)-tbsvp(jx))
      dqx = tbdsvp(jx) + (xj-jx)*(tbdsvp(jx+1)-tbdsvp(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      d = (p-erfac*qx)
      if (d .gt. 0.) then
         d    = 1.0 / d
         q    = amin1(qx * d,1.0)
         dqdt = (1.0 + erfac*q) * d * dqx
      else
         q    = 1.0
         dqdt = 0.
      endif
      return
      end
