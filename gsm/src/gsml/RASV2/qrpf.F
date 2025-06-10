      function qrpf(qrp)
!
      real qrp, qrpf, xj, c1xqrp, c2xqrp
      integer jx, nqrp
      parameter(nqrp=500001)
      real tbqrp(nqrp), tbqra(nqrp), tbqrb(nqrp)
      common/comqrp/ c1xqrp,c2xqrp,tbqrp,tbqra,tbqrb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xj   = min(max(c1xqrp+c2xqrp*qrp,1.),float(nqrp))
      jx   = min(xj,nqrp-1.)
      qrpf = tbqrp(jx)  + (xj-jx) * (tbqrp(jx+1)-tbqrp(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
