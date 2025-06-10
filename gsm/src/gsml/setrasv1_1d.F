      subroutine setrasv1_1d(lev, si, sl, del, cp, rgas ! input
     *,                 sig, sgb, prh, prj, hpk, fpk, ods, prns ! output
     *,                 lm, nstrp, rasal, krmin, krmax, afac, ufac ) ! output
!
      dimension si(lev+1), sl(lev), del(lev)
      dimension sig(lev+1), prj(lev+1), prh(lev),fpk(lev), hpk(lev)
     *,         sgb(lev),   ods(lev),   prns(lev/2), rasal(lev)
!
      rkap = rgas / cp
!
      nstrp = 1
      do l=2, lev
        if (sl(l) .gt. 0.95) nstrp = nstrp + 1
      enddo
c
      do l=1,nstrp+1
        prj(l) = (si(l)/1000.0) ** rkap
      enddo
c
      afac = 0.0
      do l=1,nstrp
        prns(l) =  (si(l)*prj(l)-si(l+1)*prj(l+1)) /
     *        ((si(l) - si(l+1)) * (rkap+1.0) )
        afac    = afac + prns(l) * del(l)
      enddo
      afac = afac / (si(1) - si(nstrp+1))
c
      lm   = lev - nstrp + 1
      lmm1 = lm - 1
      do l=1,lm
        sig(l) = si(lev-l+2)
        prj(l) = (sig(l)/1000.0) ** rkap
      enddo
      sig(lm+1) = si(1)
      prj(lm+1) = (sig(lm+1)/1000.0) ** rkap
c
      krmin = 1
      krmax = lmm1
      do l=1,lmm1
        sgb(l) = sl(lev-l+1)
        if (sgb(l) .le. 0.067) krmin = l
        if (sgb(l) .le. 0.76) krmax = l
      enddo
      sgb(lm) = 0.5 * (sig(lm) + sig(lm+1))
c
      critsig=0.95
      rasalf  = 0.30
c
      do l=1,lm
        prh(l  ) = (sig(l+1)*prj(l+1)-sig(l)*prj(l)) /
     *           ( (sig(l+1)-sig(l))*(rkap+1.0) )
        fpk(l  ) = cp * (prj(l+1)-prj(l))
        hpk(l  ) = cp * (prj(l+1)-prh(l))
        ods(l) = sig(l+1) - sig(l)
        rasal(l) = rasalf
c     
        if(sig(l).gt.critsig) rasal(l)=0.
c
      enddo
c
      do l=1,lm
        ods(l) = 1.0 / ods(l)
      enddo
      afac = 1.0 / (afac - prh(lmm1))
      ufac = 0.0
c
!
      return
      end
