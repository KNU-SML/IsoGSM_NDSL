      subroutine ddrft(
     i                  k, kd
     i,                 cp, alhl, alhf, rkap, grav, tla, alfint
     i,                 tol, qol, hol, prl, qst, hst, gam, gaf, hbl, qbl
     i,                 qrb, qrt, buy, kbl, idh, eta, rnn
     *,                 alm, wfn, train, ddft
     o,                 etd, hod, qod, evp, dof, cldfrd, wcb
     *,                 gms, gsd, ghd)
!    *,                 tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9)

!
!***********************************************************************
!******************** cumulus downdraft subroutine *********************
!****************** based on cheng and arakawa (1997)  ****** **********
!************************ subroutine ddrft  ****************************
!************************  30 april 1998  ******************************
!*******************  s. moorthi and m.j. suarez ***********************
!***********************************************************************
!***********************************************************************
!************* shrinivas.moorthi@noaa.gov (301) 763 8000(x7233) ********
!***************  max.suarez@gsfc.nasa.gov (301) 286 7373 **************
!***********************************************************************
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!===>  tol(k)     input   temperature            kelvin
!===>  qol(k)     input   specific humidity      non-dimensional

!===>  prl(k+1)   input   pressure @ edges       mb

!===>  k     input   the rise & the index of the subcloud layer
!===>  kd    input   detrainment level ( 1<= kd < k )          
!     
      implicit none
!
!  input arguments
!
      integer k, kd
      real cp, alhl, alhf, rkap, grav, alfint(k)
      integer kbl, kb1



      logical skpdd, skpup

      real hol(kd:k),   qol(kd:k),   gaf(kd:k+1)
      real hst(kd:k),   qst(kd:k),   tol(kd:k)
      real buy(kd:k+1), qrb(kd:k),   qrt(kd:k)
      real gam(kd:k+1), rnn(kd:k)
      real eta(kd:k+1), prl(kd:k+1)
!
      real    hbl,     qbl,        pris
      real    train,   wfn,        alm
!
!     temporary work space
!
      real gms(kd:k+1)
      real tx1,    tx2,  tx3, tx4
      real tx5,    tx6,  tx7, tx8, tx9
      logical unsat

      real tl, pl, ql, qs, dqs, st1,  hb, qb, tb
      real qqq, picon, piinv, del_eta
      real tem, tem1, tem2, tem3, tem4, onebg, st2
     *,    errmin, errmi2, errh, errw, erre, tem5, tem6, hbd, qbd
      integer i, l,  n, ix, kd1, ii
     *,       kp1, it, km1, ktem, kk, kk1, lm1, ll, lp1
     *,       ip1

!
      real half, one, zero, two, cmb2pa

      parameter (two=2.0,  half=0.5,  one=1.0, zero=0.0)
      parameter (cmb2pa=100.0)  !  conversion mb to pascals
!

      real        airmw, h2omw, esfac
      parameter ( airmw  = 28.97      )                                         
      parameter ( h2omw  = 18.01      )                                         
      parameter ( esfac = airmw/h2omw-one)
!
      parameter (errmin=0.0001, errmi2=0.1*errmin)
!     parameter (errmin=0.00001, errmi2=0.1*errmin)
!
      real tla, stla,  ctl2, ctl3
      real gmf, pi, vtp, onpg, rgas, ctla, vtrm, vtpexp
     *,    rpart, qrmin, aa1, bb1, cc1, dd1
     *,    wc2min, wcmin, wcbase, f2, f3, f5, gmf1, gmf5
     *,    qraf, qrbf, delx, face, angf, cmpor
!
      parameter (face=0.0, delx=10000.0, angf=0.5)
      parameter (onpg=1.0+0.5, gmf=1.0/onpg, rpart=0.0)
!     parameter (aa1=1.0, bb1=1.5, cc1=1.1, dd1=0.85, f3=cc1, f5=2.5)
      parameter (aa1=2.0, bb1=1.5, cc1=1.1, dd1=0.85, f3=cc1, f5=2.5)
!     parameter (aa1=1.0, bb1=1.0, cc1=1.0, dd1=1.0, f3=cc1,  f5=1.0)
      parameter (qrmin=1.0e-6, wc2min=0.01, gmf1=gmf/aa1, gmf5=gmf/f5)
!     parameter (qrmin=1.0e-6, wc2min=1.00, gmf1=gmf/aa1, gmf5=gmf/f5)
!
      parameter (pi=3.1415926535897931)
      integer itr, itrmu, itrmd, ktpd, itrmin
      parameter (itrmu=14, itrmd=18, itrmin=7)
!     parameter (itrmu=10, itrmd=10, itrmin=5)
      real qrp(kd:k+1), wvl(kd:k+1), al2
!
      real rnf(kd:k),   etd(kd:k+1), wcb(kd:k)
      real hod(kd:k+1), qod(kd:k+1), evp(kd:k)
      real ror(kd:k+1), stlt(kd:k)
      real ghd(kd:k),   gsd(kd:k)
      real rnt,        rnb
      real errq,       cldfrd,     rntp
      integer idn, idh
      real elm(k)
!     real em(k*k), elm(k)
      real edz, ddz, ce, qhs, fac, facg, asin, rsum1, rsum2, rsum3
      logical ddft, updret, ddlgk
!
      real aa(kd:k,kd:k+1), qw(kd:k,kd:k)
      real bud(kd:k), vt(2), vrw(2), trw(2)
     *,    gqw(kd:k)
     *,    qa(3),     wa(3),    dof, dofw
!    *,    gqw(kd:k), wcb(kd:k)

!***********************************************************************

      real qrpf, vtpf
cfpp$ expand (qrpf, qrabf, vtpf)
cfpp$ noconcur r

!


      onebg  = 1.0 / grav
      kd1    = kd + 1
      kp1    = k  + 1
      km1    = k  - 1
      kb1    = kbl - 1
!
      rgas   = rkap * cp
      cmpor  = cmb2pa / rgas
!
      vtp    = 36.34*sqrt(1.2)* (0.001)**0.1364
      vtpexp = -0.3636
      piinv  = 1.0 / pi
      picon  = pi * onebg * 0.5
!
!
!     compute rain water budget of the updraft (cheng and arakawa, 1997)
!
      cldfrd = 0.0
      rntp   = 0.0
      dof    = 0.0
      errq   = 10.0
      rnb    = 0.0
      rnt    = 0.0
      tx2    = prl(kbl)
!
      tx1      = (prl(kd) + prl(kd1)) * 0.5
      ror(kd)  = cmpor*tx1 / (tol(kd)*(1.0+esfac*qol(kd)))
!     gms(kd)  = vtp * ror(kd) ** vtpexp
      gms(kd)  = vtp * vtpf(ror(kd))
!
      qrp(kd)  = qrmin
!
      tem        = tol(k) * (1.0 + esfac * qol(k))
      ror(k+1) = 0.5 * cmpor * (prl(k+1)+prl(k)) / tem
      gms(k+1) = vtp * vtpf(ror(k+1))
      qrp(k+1) = qrmin
      buy(kd)  = max(buy(kd), 0.1)
!
      do l=kd1,k
        tem = 0.5 * (tol(l)+tol(l-1))
     *      * (1.0 + 0.5*esfac * (qol(l)+qol(l-1)))
        ror(l) = cmpor * prl(l) / tem
!       gms(l) = vtp * ror(l) ** vtpexp
        gms(l) = vtp * vtpf(ror(l))
        qrp(l) = qrmin
        buy(l) = max(buy(l), 0.1)
      enddo
!

!
      call angrad(tx1, alm, stla, ctl2, al2, pi, tla, tx2, wfn, tx3)
!
!    following ucla approach for rain profile
!
      f2      = 2.0*bb1*onebg/(pi*0.2)
      wcmin   = sqrt(wc2min)
      wcbase  = wcmin
!

      stla = f2     * stla * al2
      ctl2 = dd1    * ctl2
      ctl3 = 0.1364 * ctl2
!
      do l=kd,k
        rnf(l)   = 0.0
        wvl(l)   = 0.0
        stlt(l)  = 0.0
        gqw(l)   = 0.0
        do n=kd,k
          qw(n,l) = 0.0
        enddo
      enddo
!
!-----qw(n,l) = d(w(n)*w(n))/dqr(l)
!
      kk = kbl
      wvl(kk)    = wcbase
      qw(kd,kd)  = -qrb(kd)  * gmf1
      ghd(kd)    = eta(kd)   * eta(kd)
      gqw(kd)    = qw(kd,kd) * ghd(kd)
      gsd(kd)    = 1.0 / ghd(kd)
!
      gqw(kk)    = -  qrb(kk-1) * (gmf1+gmf1)
!
      wvl(kk)    = wcbase
      wcb(kk)    = wcbase * wcbase
      tx1        = wcb(kk)
      gsd(kk)    = 1.0
      ghd(kk)    = 1.0
!
      tem        = gmf1 + gmf1
      do l=kb1,kd1,-1
         ghd(l)  = eta(l) * eta(l)
         gsd(l)  = 1.0 / ghd(l)
         gqw(l)  = - ghd(l) * (qrb(l-1)+qrt(l)) * tem
         qw(l,l) = - qrt(l) * tem
!
         tx1     = tx1 + buy(l) * tem
         wcb(l)  = tx1 * gsd(l)
      enddo
!
      tem1        = (qrb(kd) + qrt(kd1) + qrt(kd1)) * gmf1
      gqw(kd1)    = - ghd(kd1) * tem1
!     qw(l,kd1)   = - qrt(kd1) * tem
      qw(kd1,kd1) = - qrt(kd1) * tem
      wcb(kd)     = (tx1 + buy(kd)*tem) * gsd(kd)
!
      do l=kd1,kbl
        do n=kd,l-1
           qw(n,l) = gqw(l) * gsd(n)
        enddo
      enddo
      qw(kbl,kbl) = 0.0
!
      wvl(kbl)    = wcbase
      stlt(kbl)   = 1.0 / wcbase
!
      do l=kd,k+1
        do n=kd,k
          aa(n,l) = 0.0
        enddo
      enddo
!
      skpup = .false.
!
      do itr=1,itrmu               ! rain profile iteration starts!
        if (.not. skpup) then
!
!-----calculating the vertical velocity
!
          tx1      = 0.0
          ghd(kbl) = 1.0 / qrp(kbl)
          do l=kb1,kd,-1
            tx1     = tx1    + qrp(l+1) * gqw(l+1)
            st1     = wcb(l) + qw(l,l)  * qrp(l)
     *                       + tx1      * gsd(l)
            wvl(l)  = sqrt(max(st1,wc2min))
            stlt(l) = 1.0 / wvl(l)
            ghd(l)  = 1.0 / qrp(l)
          enddo
!

!
!-----calculating trw, vrw and of
!
!         vt(1)   = gms(kd) * qrp(kd)**0.1364
          vt(1)   = gms(kd) * qrpf(qrp(kd))
          trw(1)  = eta(kd) * qrp(kd) * stlt(kd)
          tx6     = trw(1) * vt(1)
          vrw(1)  = f3*wvl(kd) - ctl2*vt(1)
          bud(kd) = stla * tx6 * qrb(kd) * 0.5
          rnf(kd) = bud(kd)
          dof     = 1.1364 * bud(kd) * ghd(kd)
          dofw    = -bud(kd) * stlt(kd)
!
          rnt     = trw(1) * vrw(1)
          tx2     = 0.0
          tx4     = 0.0
          rnb     = rnt
          tx1     = 0.5
          tx8     = 0.0
!
          if (rnt .ge. 0.0) then
            tx3 = (rnt-ctl3*tx6) * ghd(kd)
            tx5 = ctl2 * tx6 * stlt(kd)
          else
            tx3 = 0.0
            tx5 = 0.0
            rnt = 0.0
            rnb = 0.0
          endif
!
          do l=kd1,kb1
            ktem    = max(l-2, kd)
            ll      = l - 1
!
!           vt(2)   = gms(l) * qrp(l)**0.1364
            vt(2)   = gms(l) * qrpf(qrp(l))
            trw(2)  = eta(l) * qrp(l) * stlt(l)
            vrw(2)  = f3*wvl(l) - ctl2*vt(2)
            qqq     = stla * trw(2) * vt(2)
            st1     = tx1  * qrb(ll)
            bud(l)  = qqq * (st1 + qrt(l))
!
            qa(2)   = dof
            wa(2)   = dofw
            dof     = 1.1364 * bud(l) * ghd(l)
            dofw    = -bud(l) * stlt(l)
!
            rnf(ll) = rnf(ll) + qqq * st1
            rnf(l)  =           qqq * qrt(l)
!
            tem3    = vrw(1) + vrw(2)
            tem4    = trw(1) + trw(2)
!
            tx6     = .25 * tem3 * tem4
            tem4    = tem4 * ctl3
!
!-----by qr above
!
!           tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*tx7
            tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*ghd(ll)
            st1     = .25*(trw(1)*(ctl2*vt(1)-vrw(2))
     *                  * stlt(ll) + f3*trw(2))
!-----by qr below
            tem2    = .25*(trw(2)*tem3 - tem4*vt(2))*ghd(l)
            st2     = .25*(trw(2)*(ctl2*vt(2)-vrw(1))
     *                 * stlt(l)  + f3*trw(1))
!
!      from top to  the kbl-2 layer
!
            qa(1)   = tx2
            qa(2)   = qa(2) + tx3 - tem1
            qa(3)   = -tem2
!
            wa(1)   = tx4
            wa(2)   = wa(2) + tx5 - st1
            wa(3)   = -st2
!
            tx2     = tem1
            tx3     = tem2
            tx4     = st1
            tx5     = st2
!
            vt(1)   = vt(2)
            trw(1)  = trw(2)
            vrw(1)  = vrw(2)
!
            if (wvl(ktem) .eq. wcmin) wa(1) = 0.0
            if (wvl(ll)   .eq. wcmin) wa(2) = 0.0
            if (wvl(l)    .eq. wcmin) wa(3) = 0.0
            do n=ktem,kbl
              aa(ll,n) = (wa(1)*qw(ktem,n) * stlt(ktem)
     *                 +  wa(2)*qw(ll,n)   * stlt(ll)
     *                 +  wa(3)*qw(l,n)    * stlt(l) ) * 0.5
            enddo
            aa(ll,ktem) = aa(ll,ktem) + qa(1)
            aa(ll,ll)   = aa(ll,ll)   + qa(2)
            aa(ll,l)    = aa(ll,l)    + qa(3)
            bud(ll)     = (tx8 + rnn(ll)) * 0.5
     *                    - rnb + tx6 - bud(ll)
            aa(ll,kbl+1) = bud(ll)
            rnb = tx6
            tx1 = 1.0
            tx8 = rnn(ll)
          enddo
          l  = kbl
          ll = l - 1
!         vt(2)   = gms(l) * qrp(l)**0.1364
          vt(2)   = gms(l) * qrpf(qrp(l))
          trw(2)  = eta(l) * qrp(l) * stlt(l)
          vrw(2)  = f3*wvl(l) - ctl2*vt(2)
          st1     = stla * trw(2) * vt(2) * qrb(ll)
          bud(l)  = st1

          qa(2)   = dof
          wa(2)   = dofw
          dof     = 1.1364 * bud(l) * ghd(l)
          dofw    = -bud(l) * stlt(l)
!
          rnf(ll) = rnf(ll) + st1
!
          tem3    = vrw(1) + vrw(2)
          tem4    = trw(1) + trw(2)
!
          tx6     = .25 * tem3 * tem4
          tem4    = tem4 * ctl3
!
!-----by qr above
!
          tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*ghd(ll)
          st1     = .25*(trw(1)*(ctl2*vt(1)-vrw(2))
     *                * stlt(ll) + f3*trw(2))
!-----by qr below
          tem2    = .25*(trw(2)*tem3 - tem4*vt(2))*ghd(l)
          st2     = .25*(trw(2)*(ctl2*vt(2)-vrw(1))
     *                 * stlt(l)  + f3*trw(1))
!
!      for the layer next to the top of the boundary layer
!
          qa(1)   = tx2
          qa(2)   = qa(2) + tx3 - tem1
          qa(3)   = -tem2
!
          wa(1)   = tx4
          wa(2)   = wa(2) + tx5 - st1
          wa(3)   = -st2
!
          tx2     = tem1
          tx3     = tem2
          tx4     = st1
          tx5     = st2
!
          idn     = max(l-2, kd)
!
          if (wvl(idn) .eq. wcmin) wa(1) = 0.0
          if (wvl(ll)  .eq. wcmin) wa(2) = 0.0
          if (wvl(l)   .eq. wcmin) wa(3) = 0.0
!
          kk = idn
          do n=kk,l
            aa(ll,n) = (wa(1)*qw(kk,n) * stlt(kk)
     *               +  wa(2)*qw(ll,n) * stlt(ll)
     *               +  wa(3)*qw(l,n)  * stlt(l) ) * 0.5

          enddo
!
          aa(ll,idn) = aa(ll,idn) + qa(1)
          aa(ll,ll)  = aa(ll,ll)  + qa(2)
          aa(ll,l)   = aa(ll,l)   + qa(3)
          bud(ll)    = (tx8+rnn(ll)) * 0.5 - rnb + tx6 - bud(ll)
!
          aa(ll,l+1) = bud(ll)
!
          rnb        = trw(2) * vrw(2)
!
!
!      for the top of the boundary layer
!
          if (rnb .lt. 0.0) then
             kk    = kbl
             tem   = vt(2) * trw(2)
             qa(2) = (rnb - ctl3*tem) * ghd(kk)
             wa(2) = ctl2 * tem * stlt(kk)
          else
             rnb   = 0.0
             qa(2) = 0.0
             wa(2) = 0.0
          endif
!
          qa(1) = tx2
          qa(2) = dof + tx3 - qa(2)
          qa(3) = 0.0
!
          wa(1) = tx4
          wa(2) = dofw + tx5 - wa(2)
          wa(3) = 0.0
!
          kk = kbl
          if (wvl(kk-1) .eq. wcmin) wa(1) = 0.0
          if (wvl(kk)   .eq. wcmin) wa(2) = 0.0
!
          do ii=1,2
             n = kk + ii - 2
             aa(kk,n) = (wa(1)*qw(kk-1,n) * stlt(kk-1)  
     *                +  wa(2)*qw(kk,n)   * stlt(kk)) * 0.5
          enddo
          fac = 0.5
          ll  = kbl
          l   = ll + 1
          lm1 = ll - 1
          aa(ll,lm1)  = aa(ll,lm1) + qa(1)
          aa(ll,ll)   = aa(ll,ll)  + qa(2)
          bud(ll)     = 0.5*rnn(lm1) - tx6 + rnb - bud(ll)
          aa(ll,ll+1) = bud(ll)
!

!
!-----solving the budget equations for dqr
!
          do l=kd1,kbl
            lm1  = l - 1
            unsat = abs(aa(lm1,lm1)) .lt. abs(aa(l,lm1))
            do  n=lm1,kbl+1
               if (unsat) then
                  tx1       = aa(lm1,n)
                  aa(lm1,n) = aa(l,n)
                  aa(l,n)   = tx1
               endif
            enddo
            tx1 = aa(l,lm1) / aa(lm1,lm1)
            do  n=l,kbl+1
               aa(l,n) = aa(l,n) - tx1 * aa(lm1,n)
            enddo
          enddo     
!
!-----back substitution and check if the solution converges
!
          kk = kbl
          kk1 = kk + 1
          aa(kk,kk1) = aa(kk,kk1) / aa(kk,kk)      !   qr correction !
          tx2        = abs(aa(kk,kk1)) * ghd(kk)   !   error measure !
!
          kk = kbl + 1
          do l=kb1,kd,-1
             lp1   = l + 1
             tx1  = 0.0
             do n=lp1,kbl
               tx1  = tx1 + aa(l,n) * aa(n,kk)
             enddo
             aa(l,kk) = (aa(l,kk) - tx1) / aa(l,l)       ! qr correction !
             tx2      = max(tx2, abs(aa(l,kk))*ghd(l))   ! error measure !
          enddo
!
          do l=kd,kbl
             qrp(l) = max(qrp(l)+aa(l,kbl+1), qrmin)
          enddo
!
          if (itr .lt. itrmin) then
             tem = abs(errq-tx2) 
             if (tem .ge. errmi2 .and. tx2 .ge. errmin) then 
               errq  = tx2                              ! further iteration !
             else 
               skpup = .true.                           ! converges      !
               errq  = 0.0                              ! rain profile exists!
             endif 
          else
             tem = errq - tx2
             if (tem .lt. zero .and. errq .gt. 0.5) then
               skpup = .true.                           ! no convergence !
               errq = 10.0                              ! no rain profile!
             elseif (abs(tem).lt.errmi2 .or. tx2.lt.errmin) then
               skpup = .true.                           ! converges      !
               errq = 0.0                               ! rain profile exists!
             else
               errq = tx2                               ! further iteration !
             endif
          endif
!

        endif                                           ! skpup  endif!
!
      enddo                                          ! end of the itr loop!!
!

!
      if (errq .lt. 0.1) then
        ddft = .true.
        rnb  = - rnb
      else
        ddft = .false.
      endif
!
!     caution !! below is an adjustment to rain flux to maintain
!                conservation of precip!
!
      if (ddft) then
        tx1 = 0.0
        do l=kd,kb1
          tx1 = tx1 + rnf(l)
        enddo
        tx1 = train / (tx1+rnt+rnb)
        if (abs(tx1-1.0) .lt. 0.2) then
           rnt = max(rnt*tx1, 0.0)
           rnb = rnb * tx1
        else
           ddft = .false.
        endif
      endif
!
      dof = 0.0
      if (.not. ddft) return     ! rain profile did not converge!
!

      do l=kd,kb1
         rnf(l) = rnf(l) * tx1

      enddo
!
!     adjustment is over
!

!

!
!     downdraft
!
      do l=kd,k
        wcb(l) = 0.0
      enddo
!
      skpdd = .not. ddft
!
      errq  = 10.0
      if (.not. skpdd) then
!
!     calculate downdraft properties
!

        kk = max(kb1,kd1)
        do l=kk,k
          stlt(l) = stlt(l-1)
        enddo
        tem1 = 1.0 / bb1
        do l=kd,k
          evp(l) = 0.0

          if (l .lt. kbl) then
            tem     = stla * tem1
            stlt(l) = eta(l) * stlt(l) * tem / ror(l)
          else
            stlt(l) = 0.0
          endif
        enddo

!
        idn      = 99
        etd(kd)  = 0.0
        wvl(kd)  = 0.0
        hod(kd)  = hol(kd)
        qod(kd)  = qol(kd)
        tx1      = min(stlt(kd)*qrb(kd)*1.0, 1.0)    ! sigma at the top
        qrp(kd)  = 0.0
        rntp     = 0.0
        evp(kd)  = 0.0
        buy(kd)  = 0.0
        tx5      = tx1
!
!       here we assume rpart of detrained rain rnt goes to pd
!
        if (rnt .gt. 0.0) then
          rntp    = (1.0 - rpart) * rnt
          qrp(kd) = (rpart*rnt / (ror(kd)*tx1*gms(kd))) ** (1.0/1.1364)
          buy(kd) = - ror(kd) * tx1 * qrp(kd)
        endif
!
!     l-loop for the downdraft iteration from kd1 to k+1 (bottom surface)
!
!     bud(kd) = ror(kd)
      do l=kd1,k+1

          if (l .le. k) then
            st1   = 1.0 - alfint(l)
            wa(1) = alfint(l)*hol(l-1) + st1*hol(l)
            wa(2) = alfint(l)*qol(l-1) + st1*qol(l)
            wa(3) = alfint(l)*tol(l-1) + st1*tol(l)
            qa(2) = alfint(l)*hst(l-1) + st1*hst(l)
            qa(3) = alfint(l)*qst(l-1) + st1*qst(l)
          else
            wa(1) = hol(k)
            wa(2) = qol(k)
            wa(3) = tol(k)
            qa(2) = hst(k)
            qa(3) = qst(k)
          endif
!
          fac = 2.0
          if (l .eq. kd1) fac = 1.0

          facg    = fac * 0.5 * gmf5     !  12/17/97
!
          ddlgk   =  idn .eq. 99
          bud(kd) = ror(l)

          if (ddlgk) then
            tx1    = tx5
            wvl(l) = max(wvl(l-1), 0.1)

            qrp(l) = max(qrp(l-1),qrp(l))
!
!           vt(1)  = gms(l-1) * qrp(l-1) ** 0.1364
            vt(1)  = gms(l-1) * qrpf(qrp(l-1))
            rnt    = ror(l-1) * (wvl(l-1)+vt(1))*qrp(l-1)

!

            tem    = max(alm, 2.5e-4) * max(eta(l), 1.0)
            trw(1) = picon*tem*(qrb(l-1)+qrt(l-1))
            trw(2) = 1.0 / trw(1)
!
            vrw(1) = 0.5 * (gam(l-1) + gam(l))
            vrw(2) = 1.0 / (vrw(1) + vrw(1))
!
            tx4    =  (qrt(l-1)+qrb(l-1))*(onebg*fac*500.00)
!
            dofw   = 1.0 / (wa(3) * (1.0 + esfac*wa(2)))      !  1.0 / tvbar!
!
            etd(l) = etd(l-1)
            hod(l) = hod(l-1)
            qod(l) = qod(l-1)
!
            errq   = 10.0

!
            if (l .le. kbl) then
              tx3 = stlt(l-1) * qrt(l-1) * (0.5*fac)
              tx8 = stlt(l)   * qrb(l-1) * (0.5*fac)
              tx9 = tx8 + tx3
            else
              tx3 = 0.0
              tx8 = 0.0
              tx9 = 0.0
            endif
!
            tem  = wvl(l-1) + vt(1)
            if (tem .gt. 0.0) then
              tem1 = 1.0 / (tem*ror(l-1))
              tx3 = vt(1) * tem1 * ror(l-1) * tx3
              tx6 = tx1 * tem1
            else
              tx6 = 1.0
            endif
          endif
!
          if (l .eq. kd1) then
            if (rnt .gt. 0.0) then
              tem    = max(qrp(l-1),qrp(l))
              wvl(l) = tx1 * tem * qrb(l-1)*(facg*5.0)
            endif
            wvl(l) = max(0.01, wvl(l))
            trw(1) = trw(1) * 0.5
            trw(2) = trw(2) + trw(2)
          else
            if (ddlgk) evp(l-1) = evp(l-2)
          endif
!
!       no downdraft above level idh
!

          if (l .lt. idh) then

            etd(l)   = 0.0
            hod(l)   = wa(1)
            qod(l)   = wa(2)
            evp(l-1) = 0.0
            wvl(l)   = 0.0
            qrp(l)   = 0.0
            buy(l)   = 0.0
            tx5      = tx9
            errq     = 0.0
            rntp     = rntp + rnt * tx1
            rnt      = 0.0
            wcb(l-1) = 0.0
          endif
!         bud(kd) = ror(l)
!
!       iteration loop for a given level l begins
!
          do itr=1,itrmd
!
            unsat =  ddlgk .and. (errq .gt. errmin)
            if (unsat) then
!
!             vt(1)  = gms(l) * qrp(l) ** 0.1364
              vt(1)  = gms(l) * qrpf(qrp(l))
              tem    =  wvl(l) + vt(1)
!
              if (tem .gt. 0.0) then
                st1    = ror(l) * tem * qrp(l) + rnt
                if (st1 .ne. 0.0) st1 = 2.0 * evp(l-1) / st1
                tem1   = 1.0 / (tem*ror(l))
                tem2   = vt(1) * tem1 * ror(l) * tx8
              else
                tem1   = 0.0
                tem2   = tx8
                st1    = 0.0
              endif
!
              tem = ror(l)*wvl(l) - ror(l-1)*wvl(l-1)
!
              tx5   = (tx1 -st1 + tem2 + tx3) / (1.0 + tem*tem1)
              tx5   = max(tx5, 0.0)

!
              tem1   = etd(l)
              etd(l) = ror(l) * tx5 * max(wvl(l), 0.0)
!

              del_eta = etd(l) - etd(l-1)

                tem       = del_eta * trw(2)
                tem2      = max(min(tem, 1.0), -1.0)
                if (abs(tem) .gt. 1.0 .and. etd(l) .gt. 0.0 ) then
                  del_eta = tem2 * trw(1)
                  etd(l)  = etd(l-1) + del_eta
                endif
                if (wvl(l) .gt. 0.0) tx5 = etd(l) / (ror(l)*wvl(l))
                erre  = etd(l) - tem1
!

                tem1 = sqrt(max((trw(1)+del_eta)*(trw(1)-del_eta),0.0))
                edz  = (0.5 + asin(tem2)*piinv)*del_eta + tem1*piinv

              ddz   = edz - del_eta
              wcb(l-1) = etd(l) + ddz
!
              tem1  = hod(l)
              if (del_eta .gt. 0.0) then
                qqq    = 1.0 / (etd(l) + ddz)
                hod(l) = (etd(l-1)*hod(l-1) + del_eta*hol(l-1)
     *                                            + ddz*wa(1)) * qqq
                qod(l) = (etd(l-1)*qod(l-1) + del_eta*qol(l-1)
     *                                            + ddz*wa(2)) * qqq
              else
                qqq    = 1.0 / (etd(l-1) + edz)
                hod(l) = (etd(l-1)*hod(l-1) + edz*wa(1)) * qqq
                qod(l) = (etd(l-1)*qod(l-1) + edz*wa(2)) * qqq
              endif
              errh  = hod(l) - tem1
              errq  = abs(errh/hod(l))  + abs(erre/max(etd(l),0.00001))
              dof   = ddz
              vt(2) = qqq

!
              ddz  = dof
              tem4 = qod(l)
              tem1 = vrw(1)
!
              qhs  = qa(3) + 0.5 * (gaf(l-1)+gaf(l))
     *                           * (hod(l)-qa(2))
!
!                                           first iteration       !
!
              st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
              tem2 = ror(l) * qrp(l)
              call qrabf(tem2,qraf,qrbf)
              tem6 = tx5 * (1.6 + 124.9 * qraf) * qrbf * tx4
!
              ce   = tem6 * st2 / ((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!
              tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*qod(l))
              tem3   = (1.0 + tem1) * qhs * (qod(l)+ce)
              tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
              qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!

!
!                                            second iteration   !
!
              st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
              ce   = tem6 * st2 / ((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!


              tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*tem4)
              tem3   = (1.0 + tem1) * qhs * (tem4+ce)
              tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
              qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!                                              evaporation in layer l-1
!

              evp(l-1) = (qod(l)-tem4) * (etd(l)+ddz)
!                                              calculate pd (l+1/2)
              qa(1)    = tx1*rnt + rnf(l-1) - evp(l-1)
!

!
              if (etd(l) .gt. 0.0) then
                tem    = qa(1) / (etd(l)+ror(l)*tx5*vt(1))
                qrp(l) = max(tem, 0.0)
              elseif (tx5 .gt. 0.0) then
                qrp(l) = (max(0.0,qa(1)/(ror(l)*tx5*gms(l))))
     *                                          ** (1.0/1.1364)
              else
                qrp(l) = 0.0
              endif
!                                              compute buoyancy
              tem1   = wa(3)+(hod(l)-wa(1)-alhl*(qod(l)-wa(2)))
     *                                                  * (1.0/cp)
              tem1   = tem1 * (1.0 + esfac*qod(l))
              ror(l) = cmpor * prl(l) / tem1
              tem1   = tem1 * dofw
!!!           tem1   = tem1 * (1.0 + esfac*qod(l)) * dofw

              buy(l) = (tem1 - 1.0 - qrp(l)) * ror(l) * tx5
!                                              compute w (l+1/2)

              tem1   = wvl(l)
              if (etd(l) .gt. 0.0) then
              wvl(l) = vt(2) * (etd(l-1)*wvl(l-1) - facg
     *                 * (buy(l-1)*qrt(l-1)+buy(l)*qrb(l-1)))
              endif
!

!
              wvl(l) = 0.5 * (wvl(l) + tem1)

              errw   = wvl(l) - tem1
!
              errq   = errq + abs(errw/max(wvl(l), 0.00001))

!
              if (itr .ge. min(itrmin,itrmd/2)) then
                if (etd(l-1) .eq. 0.0 .and. errq .gt. 0.2) then
                  ror(l)   = bud(kd)
                  etd(l)   = 0.0
                  wvl(l)   = 0.0
                  errq     = 0.0
                  hod(l)   = wa(1)
                  qod(l)   = wa(2)
                  tx5      = tx1 + tx9

                  evp(l-1) = 0.0
                  tem      = max(tx1*rnt+rnf(l-1), 0.0)
                  qa(1)    = tem - evp(l-1)
                  qrp(l)   = (qa(1) / (ror(l)*tx5*gms(l)))
     *                                            ** (1.0/1.1364)
                  buy(l)   = - ror(l) * tx5 * qrp(l)
                  wcb(l-1) = 0.0
                endif
!
                del_eta = etd(l) - etd(l-1)
                if(del_eta .lt. 0.0 .and. errq .gt. 0.1) then
                  ror(l)   = bud(kd)
                  etd(l)   = 0.0
                  wvl(l)   = 0.0
                  tx5      = tx1 + tx9
                  cldfrd   = tx5
!
                  del_eta  = - etd(l-1)
                  edz      = 0.0
                  ddz      = -del_eta
                  wcb(l-1) = ddz

!
                  hod(l)   = hod(l-1)
                  qod(l)   = qod(l-1)

!
                  tem4     = qod(l)
                  tem1     = vrw(1)
!
                  qhs      = qa(3) + 0.5 * (gaf(l-1)+gaf(l))
     *                                   * (hod(l)-qa(2))

!
!                                           first iteration       !
!
                  st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
                  tem2 = ror(l) * qrp(l-1)
                  call qrabf(tem2,qraf,qrbf)
                  tem6 = tx5 * (1.6 + 124.9 * qraf) * qrbf * tx4
!
                  ce   = tem6*st2/((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!

                  tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*qod(l))
                  tem3   = (1.0 + tem1) * qhs * (qod(l)+ce)
                  tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
                  qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!
!                                            second iteration   !
!
                  st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
                  ce   = tem6*st2/((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!


                  tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*tem4)
                  tem3   = (1.0 + tem1) * qhs * (tem4+ce)
                  tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
                  qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))

!                                              evaporation in layer l-1
!
                  evp(l-1) = (qod(l)-tem4) * (etd(l)+ddz)

!                                               calculate pd (l+1/2)
                  rnn(l-1) = tx1*rnt + rnf(l-1) - evp(l-1)
 
                  errq = 0.0
                endif
              endif
            endif
!

          enddo                ! end of the iteration loop  for a given l!
          if (l .le. k) then
            if (etd(l-1) .eq. 0.0
     *         .and. errq .gt. errmin*10.0) then
               ror(l)   = bud(kd)
               tx5      = tx1 + tx9
               evp(l-1) = 0.0
               qrp(l)   = 0.0
               qa(1)    = tx1*rnt + rnf(l-1) - evp(l-1)
               if (qa(1) .gt. 0.0) then
                 qrp(l) = (qa(1) / (ror(l)*tx5*gms(l)))
     *                                         ** (1.0/1.1364)
               endif
               etd(l)   = 0.0
               wvl(l)   = 0.0
               st1      = 1.0 - alfint(l)
               hod(l)   = wa(1)
               qod(l)   = wa(2)

               errq     = 0.0
               buy(l)   = - ror(l) * tx5 * qrp(l)
               wcb(l-1) = 0.0
            endif
          endif

!
!
          ll = min(idn, k+1)
          if (errq .lt. 1.0 .and. l .le. ll) then
            if (etd(l-1) .gt. 0.0 .and. etd(l) .eq. 0.0) then
             idn = l
            endif
          endif


! 
!     if downdraft properties are not obtainable, (i.e.solution does
!      not converge) , no downdraft is assumed
!
           if (errq .gt. errmin*10.0 .and. idn .eq. 99)
     *                          ddft = .false.
!
!
        dof = 0.0
        if (.not. ddft) return
!

      enddo                      ! end of the l loop of downdraft !

        tx1 = 0.0

        dof = qa(1)
!

      endif                       ! skpdd endif
!

       rnn(kd) = rntp
       tx1     = evp(kd)
       tx2     = rntp + rnb

       ii = idh
       if (ii .ge. kd1+1) then
          rnn(kd)   = rnn(kd) + rnf(kd)
          tx2       = tx2 + rnf(kd)
          rnn(ii-1) = 0.0
          tx1       = evp(ii-1)
        endif
        do l=kd,k
          ii = idh

          if (l .gt. kd1 .and. l .lt. ii) then
            rnn(l-1) = rnf(l-1)
            tx2      = tx2 + rnn(l-1)

          elseif (l .ge. ii .and. l .lt. idn) then

            rnn(l)   = 0.0
            tx1      = tx1 + evp(l)

          elseif (l .gt. idn) then
            etd(l)   = 0.0
            hod(l)   = 0.0
            qod(l)   = 0.0
            evp(l-1) = 0.0
            rnn(l-1) = rnf(l-1)
            tx2      = tx2 + rnn(l-1)
          endif
        enddo
        if (k+1 .gt. idn) then
          etd(k+1) = 0.0
          hod(k+1) = 0.0
          qod(k+1) = 0.0
          evp(k)   = 0.0
          rnn(k)   = rnf(k)
        endif
!
!      for downdraft case the rain is that falls thru the bottom

        l = kbl

        rnn(l)  = rnb
        cldfrd  = tx5

        tx2     = tx2 + dof

!
!     caution !! below is an adjustment to rain flux to maintain
!                conservation of precip!

!

        if (tx1 .gt. 0.0) then
          tx1 = (train - tx2) / tx1
        else
          tx1 = 0.0
        endif

        tx5      = evp(kbl)

        evp(kbl) = evp(kbl) * tx1

        tx3      = rnn(kbl) + evp(kbl) + dof
        tx2      = rnn(kbl)
        tx4      = evp(kbl)


        do l=kd,kb1

          tx5    = tx5 + evp(l)
          evp(l) = evp(l) * tx1
          tx3    = tx3 + evp(l) + rnn(l)
          tx2    = tx2 + rnn(l)
          tx4    = tx4 + evp(l)
        enddo


!
!***********************************************************************
!***********************************************************************

      return
      end
