      subroutine setrasv2_hyb(lm, prsi, prsl, cp, rgas
     *,                 sig, sgb, prj
     *,                 krmin, krmax, kfmax
     *,                 mct, kctop)
!
      dimension prsi(lm+1),prsl(lm)
      dimension sig(lm+1), sgb(lm),   prj(lm+1)
!
      dimension pctop(mct+1), kctop(mct+1)
      data pctp/70.0/, pcbt/900.0/, pinc/5.0/ 
!     should be consistent with setrasv2.F
!
      rkap = rgas / cp
      lmm1 = lm - 1
      do l=1,lm
         sig(l) = prsi(lm-l+2)/prsi(1)
         prj(l) = (0.001*sig(l)) ** rkap
         sgb(l) = prsl(lm-l+1)/prsi(1)
      enddo
      sig(lm+1) = prsi(1)/prsi(1)
      prj(lm+1) = (0.001*sig(lm+1)) ** rkap
c
      tem  = (pcbt - pctp - (mct-1)*mct*pinc/2) / mct
      pctop(1) = pctp * 0.001
      do i=1,mct
         pctop(i+1) = 0.001 * (pctop(i)*1000.0 + tem + pinc*(i-1))
      enddo
      do l=1,lm
         if (sig(l) . le. pctop(1)) kctop(1) = l
      enddo
      tem2 = 0.0
      do i=2,mct+1
         if (kctop(i-1) .lt. lm) then
            ii       = kctop(i-1) + 1
            kctop(i) = ii
            do l=ii,lm
               if (sig(l) . le. pctop(i)) kctop(i) = l
            enddo
            if (pctop(i) .gt. sig(kctop(i))) then
               tem1 = sig(kctop(i)) - sig(kctop(i-1))
               if (tem1 .lt. tem2) kctop(i) = min(lm, kctop(i)+1)
               tem2 = sig(kctop(i)) - sig(kctop(i-1))
            endif
         else
            kctop(i) = lm + 1
         endif
      enddo
      do i=1,mct+1
         if (kctop(i) .le. lm) mctm = i
      enddo
c
      krmin = 1
      krmax = mctm - 1
      kfmax = mctm - 1
      do l=1,mctm-1
         sgc  = 0.5 * (sig(kctop(l)) + sig(kctop(l+1)))
         if (sgc .le. 0.760) krmax = l
         if (sgc .le. 0.930) kfmax = l
      enddo      
!
      return
      end
