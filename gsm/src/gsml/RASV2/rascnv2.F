      subroutine rascnv2(lm, dt
     *,           krmin, krmax, kfmax, ncrnd, rannum, pdd
     *,           mct, kctop
!    *,           mct, kctop, sgh
     *,           rgas, cp, grav, alhl, alhf
     *,           tin, qin, uin, vin, pin, kpbl
     *,           rainc, kbot, ktop, kuo, lat, cdrag
     *,           sig, prj, sgb, rasal, clw, clt, dsfc) 
!
!
      implicit none
      logical revap, cumfrc
      logical wrkfun, calkbl, crtfun, updret, lprnt
      real frac, crtmsf, tf, tcr, tcrf, max_neg_bouy
      parameter (frac=0.5, crtmsf=0.0)
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
      parameter (max_neg_bouy=0.20, revap=.true., cumfrc=.false.)
      parameter (wrkfun = .false.,  updret = .false.)
!     parameter (crtfun = .false.,  calkbl = .true.)
      parameter (crtfun = .true.,   calkbl = .true.)
!
!      input
!
      integer lm, mct, lat, kbot, ktop, kuo, kbl, kpbl
     *,       krmin, krmax, kfmax, ncrnd
      real    tin(lm), qin(lm), uin(lm), vin(lm), pin
     *,       rainc,  clw(lm), cdrag, rannum(ncrnd), clt(lm), dsfc
     *,       dt, pdd, rgas, cp, grav, alhl, alhf
!
      real sig(lm+1),  prj(lm+1), alfint(lm), sgb(lm), rasal(lm)
      integer         kctop(mct+1)
!
!     locals
!
      real  pib,  pik, rain,  uvc(lm,2), tn0(lm), qn0(lm) 
     *,     tcu(lm),   qcu(lm), pcu(lm)
     *,     flx(lm),   qii(lm), qli(lm),  prs(lm+1), psj(lm+1)
     *,     rcu(lm,2), rkap
      integer i, ia, l, ll, lmm1, lmp1
!
!     logical lprnt
!
cfpp$ noconcur r
!
!     print *,' sgb=',sgb
      lmp1 = lm + 1
!
!
      lmm1 = lm - 1
      lmp1 = lm + 1
      rkap = rgas / cp
!
      alfint(1) = 0.5
      do l=2,lm
!        alfint(l) = (sgb(l) - sig(l)) / (sgb(l) - sgb(l-1))
         alfint(l) = 0.5
      enddo
!     write(0,*)' mct=',mct,' kctop=',kctop,' sgh=',sgh
!
      ia = 1
!
!     print *,' in rascnv: lm=',lm,'lat=',lat,' ia=',ia
!     print *,' jr=',jr
!    *,' lprnt=',lprnt
!     if (lprnt) then
!        if (me .eq. 0) then
!        print *,' tin',(tin(ia,l),l=lm,1,-1)
!        print *,' qin',(qin(ia,l),l=lm,1,-1)
!     endif
!
!     assumes initial value of cloud liquid water to be zero
!     i.e. no environmental liquid water used!!!
!
      do l=1,lm
         clw(l) = 0.0
         clt(l) = 0.0
         qii(l)   = 0.0
         qli(l)   = 0.0
!
!      initialize cloudiness,heating, drying, cloudiness etc.
!
         tcu(l)   = 0.0
         qcu(l)   = 0.0
         pcu(l)   = 0.0
         flx(l)   = 0.0
         rcu(l,1) = 0.0
         rcu(l,2) = 0.0
      enddo
      pib = pin * 10.0
      pik = pib ** rkap
!
!     transfer input prognostic data in to local variable
!
      do l=1,lm
         ll = lm - l + 1
         tn0(l)   = tin(ll)
         qn0(l)   = qin(ll)
         uvc(l,1) = uin(ll)
         uvc(l,2) = vin(ll)
! 10/21/97
!        qii(l) = clw(l) 
!    *          * max(0.0, min(1.0, (tcr-tn0(l))*tcrf))
!        qli(l) = clw(l) - qii(l)
! 10/21/97
!
      enddo
      kbl  = max(min(lm, lm+1-kpbl), lm/2)
      rain = 0.0
!
      do l=1,lmp1
         prs(l) = pib * sig(l)
         psj(l) = pik * prj(l)
      enddo
!
!     write(0,*)' calling rasv2: ' lm=',lm,' dt=',dt
!    *,' lat=',lat
!     if (lprnt) then
!        print *,' qns=',qns(ia),' qn0=',qn0(ia,lm),'qin=',qin(ia,1)
!    *,' pik=',pik(ia)
!        if (me .eq. 0) then
!        print *,' tn0',(tn0(ia,l),l=1,lm)
!        print *,' qn0',(qn0(ia,l),l=1,lm),' kbl=',kbl
!     endif
!
!     print *,' alhf=',alhf,' max=',max_neg_bouy,' alfint=',alfint
!
      call rasv2(lm,  dt, pdd
     *,          ncrnd, krmin, krmax, kfmax, frac, rasal, .true.
     *,          revap, max_neg_bouy, alfint, cumfrc
!    *,          revap, max_neg_bouy, alfint, alfinq
     *,          cp,  alhl, alhf, grav, rkap, rannum
     *,          mct, kctop
!    *,          mct, kctop, sgh, sig
!
     *,          tn0, qn0, uvc, kbl
!    *,          qli,    qii
     *,          qli,    qii, clw
     *,          rain,   clt, tcu, qcu, pcu, flx, rcu
     *,          cdrag,  dsfc
     *,          prs,    psj, wrkfun, calkbl, crtfun, updret)
!
!
#ifdef DBG
!     tem1 = 0.0
!     tem2 = 0.0
!     tem3 = 0.0
!     do l=1,lm
!        cf1 = (100.0*cp)*(sig(l+1)-sig(l)) / (alhl*grav)
!        cf2 = 100.0 *(sig(l+1)-sig(l))/ grav
!        print *,' cp=',cp,' alhl=',alhl,' l=',l
!    *,' cf1=',cf1,' cf2=',cf2
!        print *,' tcu=',tcu(l)
!        print *,' qcu=',qcu(l)
!           tem1 = tem1 + tcu(l)*pib*cf1
!           tem2 = tem2 + qcu(l)*pib*cf2
!     enddo
!     tem3 = tem3 + rain
!     print *,' after rasv2: tem1=',tem1,' tem2=',tem2,' tem3=',tem3
!     if (tem1 .ne. 0.0) stop

!
!     if (lprnt) then
!        print *,' tn0',(tn0(ia,l),l=1,lm)
!        print *,' qn0',(qn0(ia,l),l=1,lm)
!        print *,' qns=',qns(ia),' qn0=',qn0(ia,lm),'qin=',qin(ia,1)
!     endif
!
#endif
!
        rainc = rain * 0.001
!
!     if (lprnt) then
!        print *,' tn0',(tn0(imax,l),l=1,lm)
!        print *,' qn0',(qn0(imax,l),l=1,lm)
!     endif
!
      do l=1,lm
         ll = lm - l + 1
         tin(ll) = tn0(l)
         qin(ll) = qn0(l)
         uin(ll) = uvc(l,1)
         vin(ll) = uvc(l,2)
!        clw(l ) = qli(l) + qii(i,l)  ! now done in rasv2!
      enddo
      kuo  = 0
!     kbot = lmp1
!     ktop = 0
      ktop = lmp1
      kbot = 0

      do l=lmm1,1,-1
        if (tcu(l) .ne. 0.0) then
           kuo = 1
        endif
      enddo
!  new test for convective clouds ! added in 08/21/96
      do l=lmm1,1,-1
        if (clw(l) .gt. 0.0) ktop = l
      enddo
      do l=1,lmm1
        if (clw(l) .gt. 0.0) kbot = l
      enddo
      ktop = lmp1 - ktop
      kbot = lmp1 - kbot
!
!     if (lprnt) then
!        print *,' tin',(tin(ia,l),l=lm,1,-1)
!        print *,' qin',(qin(ia,l),l=lm,1,-1)
!     endif
!
      return
      end
