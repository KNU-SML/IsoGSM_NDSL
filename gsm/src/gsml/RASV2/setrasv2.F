      subroutine setrasv2(lm, si, sl, cp, rgas, dt, nsphys, fh
     *,                 sig, sgb, prj
     *,                 rasal, krmin, krmax, kfmax, ncrnd, rannum
     *,                 mct, kctop, dtf, dsfc, len)
!
!     this version of setras assumes same number of levels in the vertical
!     but levels are reverse order
!
!
      dimension si(lm+1),  sl(lm),    dsfc(len)
      dimension sig(lm+1), sgb(lm),   prj(lm+1)
     *,         rasal(lm), rannum(200,nsphys)
!
      dimension  kctop(mct+1)
      data pctp/70.0/, pcbt/900.0/, pinc/5.0/
!
      logical first
!
cfpp$ noconcur r
!
      call sigset(lm, si, sl, cp, rgas, sig, sgb, prj)
      call rasset(lm, dt, nsphys, fh,  sig
     *,           rasal, krmin, krmax, kfmax, ncrnd, rannum
     *,           mct, pctp, pcbt, kctop, dtf, dsfc, len)
!
      return
      end
