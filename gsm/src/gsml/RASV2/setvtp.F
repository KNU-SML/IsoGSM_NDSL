      subroutine setvtp
      parameter(nvtp=10001, vtpexp=-0.3636)
      real tbvtp(nvtp)
      common/comvtp/ c1xvtp,c2xvtp,tbvtp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfpp$ noconcur r
      xmin=0.05
      xmax=1.5
      xinc=(xmax-xmin)/(nvtp-1)
      c1xvtp=1.-xmin/xinc
      c2xvtp=1./xinc
      do jx=1,nvtp
        x         = xmin + (jx-1)*xinc
        tbvtp(jx) =        x ** vtpexp
      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
