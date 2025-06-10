      subroutine qrabf(qrp,qraf,qrbf)
!
      real qrp, qraf, qrbf, xj, c1xqrp, c2xqrp
      integer jx, nqrp
      parameter(nqrp=500001)
      real tbqrp(nqrp), tbqra(nqrp), tbqrb(nqrp)
      common/comqrp/ c1xqrp,c2xqrp,tbqrp,tbqra,tbqrb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xj   = min(max(c1xqrp+c2xqrp*qrp,1.),float(nqrp))
      jx   = min(xj,nqrp-1.)
      xj   = xj - jx
      qraf = tbqra(jx)  + xj * (tbqra(jx+1)-tbqra(jx))
      qrbf = tbqrb(jx)  + xj * (tbqrb(jx+1)-tbqrb(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
