      subroutine sigset(lm, si, sl, cp, rgas, sig, sgb, prj)
!
!     this version of setras assumes same number of levels in the vertical
!     but levels are reverse order
!
!
      dimension si(lm+1),  sl(lm)
      dimension sig(lm+1), sgb(lm),   prj(lm+1)
!
      logical first
      data first/.true./
      save first
!
cfpp$ noconcur r
!
      if (first) then
         rkap = rgas / cp
         lmm1 = lm - 1
         do l=1,lm
            sig(l) = si(lm-l+2)
            prj(l) = (0.001*sig(l)) ** rkap
            sgb(l) = sl(lm-l+1)
         enddo
         sig(lm+1) = si(1)
         prj(lm+1) = (0.001*sig(lm+1)) ** rkap
!
      endif
!
      return
      end
