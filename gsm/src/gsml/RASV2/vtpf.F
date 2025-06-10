      function vtpf(ror)
!
      real ror, vtpf, xj, c1xvtp, c2xvtp
      integer jx, nvtp
      parameter(nvtp=10001)
      real tbvtp(nvtp)
      common/comvtp/ c1xvtp,c2xvtp,tbvtp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xj   = min(max(c1xvtp+c2xvtp*ror,1.),float(nvtp))
      jx   = min(xj,nvtp-1.)
      vtpf = tbvtp(jx)  + (xj-jx) * (tbvtp(jx+1)-tbvtp(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
