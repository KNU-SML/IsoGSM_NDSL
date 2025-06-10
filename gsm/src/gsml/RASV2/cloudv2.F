      subroutine cloudv2(

     i                  k, kd, m
     i,                 rasalf, fracbl, max_neg_bouy
     i,                 alfint, alfinq, rhfacl, rhfacs
     i,                 cp, alhl, alhf, rkap, grav

     i,                 toi, qoi, roi, prs, prj
     i,                 qli, qii, kpbl, dsfc
     i,                 cd
     i,                 idiag

     u,                 tcu, qcu, rcu, pcu, flx
!    u,                 tcd, qcd
     u,                 cup, revap, dt
     *,                 wfnc, wrkfun, calkbl, crtfun, tla, dndrft, pdd)

!
!***********************************************************************
!******************** relaxed  arakawa-schubert ************************
!****************** plug compatible scalar version *********************
!************************ subroutine cloud  ****************************
!************************  23 june 1999     ****************************
!**************************  version 2.0  ******************************
!*******************  s. moorthi and m.j. suarez ***********************
!***********************************************************************
!*reference:
!-----------
!     noaa technical report nws/ncep 99-01:
!     documentation of version 2 of relaxed-arakawa-schubert
!     cumulus parameterization with convective downdrafts, june 1999.
!
!***********************************************************************
!
!===>    updates cloud tendencies due to a single cloud
!===>    detraining at level kd.
!
!***********************************************************************
!************* shrinivas.moorthi@noaa.gov (301) 763 8000(x7233) ********
!***************  max.suarez@gsfc.nasa.gov (301) 286 7373 **************
!***********************************************************************
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!===>  toi(k)     inout   temperature            kelvin
!===>  qoi(k)     inout   specific humidity      non-dimensional
!===>  roi(k,m)   inout   tracer                 arbitrary
!===>  qli(k)     inout   liquid water           non-dimensional
!===>  qii(k)     inout   ice                    non-dimensional

!===>  prs(k+1)   input   pressure @ edges       mb
!===>  prj(k+1)   input   (p/p0)^kappa  @ edges  non-dimensional

!===>  k      input   the rise & the index of the subcloud layer
!===>  kd     input   detrainment level ( 1<= kd < k )          
!===>  m      input   number of tracers. may be zero.
!===>  dndrft input   logical .true. or .false.
!===>  pdd    input   pressure level above which downdrfat can exist hpa
!
!===>  tcu(k  )   update  temperature tendency       deg
!===>  qcu(k  )   update  water vapor tendency       (g/g)
!===>  rcu(k,m)   update  tracer tendencies          nd
!===>  pcu(k-1)   update  precip @ base of layer     kg/m^2
!===>  flx(k  )   update  mass flux @ top of layer   kg/m^2
!===>  cup        update  precipitation at the surface kg/m^2
!
      implicit none
!
      real almin1, almin2, almax
      parameter (almin1=0.00e-6, almin2=0.00e-5, almax=1.0e-1)
!     parameter (almin1=0.10e-4, almin2=0.15e-4, almax=1.0e-1)

!  input arguments

      logical revap, dndrft, wrkfun, calkbl, crtfun, calcup
      integer k, kd, m

      integer idiag

      real toi(k  ),  qoi(k ),  prs(k+1)
      real qli(k),    qii(k)
      real prj(k+1),  roi(k,m)
      real cd,        ufn,     dsfc
      integer kpbl,   kbl,     kb1

      real rasalf, fracbl, max_neg_bouy, alfint(k), rhfacl, rhfacs
      real cp, alhl, alhf, rkap, grav, alfinq(k),   pdd

 
!  update arguments

      real tcu(k), qcu(k), rcu(k,m)
      real tcd(k), qcd(k)
      real pcu(k), flx(k)
      real cup

!  temporary work space

      real hol(kd:k),  qol(kd:k),   gaf(kd:k+1)
      real hst(kd:k),  qst(kd:k),   tol(kd:k)
      real gmh(kd:k),  gms(kd:k+1), gam(kd:k+1)
      real akt(kd:k),  akc(kd:k),   bkc(kd:k)
      real ltl(kd:k),  rnn(kd:k),   fco(kd:k)
      real prh(kd:k),  pri(kd:k)
      real qil(kd:k),  qll(kd:k)
      real dlb(kd:k+1),dlt(kd:k+1), eta(kd:k+1)
      real prl(kd:k+1)

      real alm,   det,    hcc,  clp
      real hsu,   hsd,    qtl,  qtv
      real akm,   wfn,    hos,  qos
      real amb,   tx1,    tx2,  tx3
      real tx4,   tx5,    qis,  qls
      real hbl,   qbl,    rbl(m)
      real qlb,   qib,    pris
      real wfnc,  tx6,    acr
      real tx7,   tx8,    tx9

      integer ia,  i1,  i2, id1, id2
      integer ib,  i3

      logical unsat

      logical lowest, skpdd

      real tl, pl, ql, qs, dqs, st1, sgn, c0, tau, qtvp, hb, qb, tb
      real hccp, ds, dh, ambmax, x00, epp, qtlp, qqq, c0i
      real dpi, dphib, dphit, del_eta, detp, qudfac, testmb, hcrit
      real tem, tem1, tem2, tem3, tem4, onebg, rkpp1i, st2
     *,    errmin, errmi2, errh, errw, erre, tem5, tem6, hbd, qbd
      integer i, l,  lend1, lend2, n,  kd1, ii
     *,       kp1, it, km1, ktem, kk, kk1, lm1, ll, lp1



!     reevaporation

       real clfrac, dt

       real actevap,arearat,deltaq,mass,massinv,potevap
       real teq,qsteq,dqdt,qeq
       real elocp,gravcon,afc, rknob, elfocp
!
      real half, one, zero, two, cmb2pa, rhmax, quad_lam, rhram

      parameter (two=2.0,  half=0.5,  one=1.0, zero=0.0)

      parameter (cmb2pa=100.0)  !  conversion mb to pascals
      parameter (rhmax=1.0   )  !  max relative humidity
      parameter (quad_lam=1.0)  !  mask for quadratic lambda
      parameter (rhram=0.15)    !  pbl relative humidity ramp
      parameter (hcrit=3000.0)  !  critical moist static energy
!

      real        airmw, h2omw, esfac
      parameter ( airmw  = 28.97      )                                         
      parameter ( h2omw  = 18.01      )                                         
      parameter ( esfac = airmw/h2omw-one)
!
      integer kblmx
      parameter (c0i=2.0e-3)
      parameter (c0=2.0e-3, kblmx=10, errmin=0.0001, errmi2=0.1*errmin)
!     parameter (c0=0.5e-3, kblmx=7, errmin=0.00001, errmi2=0.1*errmin)
      real tf, tcr, tcrf
      parameter (tf=230.16, tcr=260.16, tcrf=1.0/(tcr-tf))
!     parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!     temporary workspace and parameters needed for downdraft
!
      real tla, gmf
!
!     parameter (face=0.0, delx=10000.0, angf=0.5)
!

      real buy(kd:k+1), qrb(kd:k),   qrt(kd:k)
      real etd(kd:k+1), hod(kd:k+1), qod(kd:k+1)
      real ghd(kd:k),   gsd(kd:k),   evp(kd:k)
      real etz(kd:k)
      real train, dof, cldfrd
      integer idh
      real fac, rsum1, rsum2, rsum3
      logical ddft, updret

!

!***********************************************************************

cfpp$ expand (qsatcn, crtwrk)
cfpp$ noconcur r

      tcd = 0.0
      qcd = 0.0
!

      qudfac = quad_lam*half
      onebg  = 1.0 / grav
      rkpp1i = one / (one+rkap)
      kd1    = kd + 1
      kp1    = k  + 1
      km1    = k  - 1
!
      do l=kd, k
        rnn(l) = 0.0
      enddo
!
      elocp  = alhl / cp
      elfocp = (alhl+alhf) / cp
!

!
      dlt(kp1) = prj(kp1)
      cldfrd   = 0.0
      dof      = 0.0
      prl(kp1) = prs(kp1)
      do l=kd,k
        dlt(l) = prj(l)
        tol(l) = toi(l)
        qol(l) = qoi(l)
        prl(l) = prs(l)
        buy(l) = 0.0
      enddo
!
      do l=kd, k
 
        dpi = one / (prl(l+1) - prl(l))
        pri(l) = (grav/cmb2pa) * dpi
 
        prh(l) = (dlt(l+1)*prl(l+1) - dlt(l)*prl(l))*dpi * rkpp1i
 
        pl     = (prl(l+1) + prl(l))*half
!       tl     = pol(l)*prh(l)
        tl     = tol(l)
 
        akt(l) = (prl(l+1) - pl) * dpi
!       akt(l) = half
 
        call qsatcn(tl, pl, qs, dqs)
 
        qst(l) = qs
        gam(l) = dqs * elocp
        gaf(l) = (one/alhl)*(gam(l)/(one+gam(l)))
 
        ql     = amax1(amin1(qs*rhmax,qol(l)), 1.0e-10)
        qol(l) = ql
 
        tem    = cp * tl
        ltl(l) = tem*(one+gam(l))/(one+esfac*(qst(l)+tl*dqs))
        eta(l) = one / (ltl(l) * (1.0 + esfac*ql))

        hol(l) = tem + ql * alhl
        hst(l) = tem + qs * alhl

      enddo
!
      eta(k+1) = zero
      gms(k)   = zero
!
      akt(kd)  = half
      gms(kd)  = zero
!
      clp      = zero
!
      gam(k+1) = gam(k)
      gaf(k+1) = gaf(k)
!
      do l=k,kd1,-1
        tem    = cp * tol(l) / prh(l)
        tem1   = (one + esfac*qol(l)) * tem
 
        dphib  = (prj(l+1) - prh(l)) * tem1
        dphit  = (prh(l  ) - prj(l)) * tem1
!
        tem    = tem * eta(l)
        dlb(l) = (prj(l+1) - prh(l)) * tem
        dlt(l) = (prh(l  ) - prj(l)) * tem
!
        qrb(l) = dphib
        qrt(l) = dphit
!
        eta(l) = eta(l+1) + dphib

        hol(l) = hol(l) + eta(l)
        hst(l) = hst(l) + eta(l)

        eta(l) = eta(l) + dphit
      enddo
!
!     for the cloud top layer
!
      l = kd
      tem    = cp * tol(l) * (prj(l+1) - prh(l)) / prh(l)
      dphib  = tem * (one + esfac*qol(l))
!
      dlb(l) = tem * eta(l)
!
      qrb(l) = dphib
      qrt(l) = dphib
!
      eta(l) = eta(l+1) + dphib

      hol(l) = hol(l) + eta(l)
      hst(l) = hst(l) + eta(l)
!
!     modification of the environmental moist static energy if ice present
!
      do l=kd,k
        hol(l) = hol(l) - alhf * qii(l)
      enddo
!
!     to determine kbl internally -- if kbl is defined externally
!     the following two loop should be skipped
!
      if (calkbl) then
         kbl   = k
         tx1   = zero
         unsat = .false.
         ktem = max(kd, k-kblmx-2)
         do l=km1,ktem,-1
           tem = hol(k) - hol(l)
           tx3 = (hol(l) - hol(l+1)) / (prl(l+2) - prl(l))

           if (tx3 .lt. tx1 .and. tem .lt. hcrit) then
             tx1   = tx3
             kbl   = l+1
             unsat = .true.
           elseif (unsat .and.
     *           ( ((kbl .lt. k-1) .and. tx3 .gt. 0.5*tx1)
     *              .or. tem .gt. hcrit) ) then
             tx1 = -1.0e20
           endif
         enddo
!
         do l=ktem,k
           if (hol(k) .gt. hst(l)) kbl = l
         enddo
         kbl = min(k, max(kbl+1, kd-1))
!
         kbl  = max(kbl,k-kblmx)

         kpbl = kbl
      else
         kbl  = kpbl
      endif
!

      kbl      = max(kbl,kd)
      kb1      = kbl - 1
      tx1      = eta(kbl)
!
      pris     = one / (prl(k+1)-prl(kbl))
!
      gms(kbl) = 0.0
      gmh(kbl) = 0.0
      rnn(kbl) = 0.0

!
      do l=k,kd,-1
        if (l .ge. kbl) then
          eta(l) = (prl(k+1)-prl(l)) * pris
        else
          rnn(l) = (eta(l) - tx1) * onebg
          eta(l) =  rnn(l) - rnn(l+1)
          gmh(l) =  rnn(l) * rnn(l) * qudfac
          gms(l) =  gmh(l) - gmh(l+1)
        endif
      enddo
!
      hbl = hol(k) * eta(k)
      qbl = qol(k) * eta(k)
      qlb = qli(k) * eta(k)
      qib = qii(k) * eta(k)
      tx1 = qst(k) * eta(k)
!
      do l=km1,kbl,-1
         tem = eta(l) - eta(l+1)
         hbl = hbl + hol(l) * tem
         qbl = qbl + qol(l) * tem
         qlb = qlb + qli(l) * tem
         qib = qib + qii(l) * tem
         tx1 = tx1 + qst(l) * tem
      enddo
!                                   find min value of hol in tx2
      tx2 = hol(kd)
      idh = kd + 1
      do l=kd1,kb1
        if (hol(l) .lt. tx2) then
           tx2 = hol(l)
           idh = l             ! level of minimum moist static energy!
        endif
      enddo
!     idh = 1
      idh = max(kd+1, idh)
!
      tem1 = hbl - hol(kd)
      tem  = hbl - hst(kd1)
     *             - ltl(kd1) *( esfac *(qol(kd1)-qst(kd1)))
      lowest = kd .eq. kb1
!

      tx1   = qbl / tx1
      unsat = (tem .gt. zero .or. (lowest .and. tem1 .ge. zero))
     *         .and. (tx1 .gt. rhfacs-rhram)
     *         .and. (kbl .gt. kd)

!
!===>  if no sounding meets first condition, return

      if (.not. unsat) return
!
!
      do n=1,m
        rbl(n) = roi(k,n) * eta(k)
        do l=km1,kbl,-1
          rbl(n) = rbl(n) + roi(l,n)*(eta(l)-eta(l+1))
        enddo
      enddo
!                   linear case!
!     tem = 1.0 / rhram
!     qtl    = max(zero, min(1.0, (tx1-rhfacs+rhram)*tem))
!
!  commented on 09/23/98 to test linear case!
      tem1   = tx1 - rhfacs
      qtl    = max(zero, min(one, exp(20.0*tem1) ))
!

      tx4    = 0.0
      tx5    = 0.0
!
      tx3      = qst(kbl) - gaf(kbl) * hst(kbl)
      qil(kbl) = max(0.0, min(1.0, (tcr-tol(kbl))*tcrf))
!
      do l=kb1,kd1,-1
        tem      = qst(l) - gaf(l) * hst(l)
        tem1     = (tx3 + tem) * 0.5
        st2      = (gaf(l)+gaf(l+1)) * 0.5
!
        fco(l+1) =            tem1 + st2 * hbl
        rnn(l+1) = rnn(l+1) * tem1 + st2 * tx4
        gmh(l+1) = gmh(l+1) * tem1 + st2 * tx5
!
        tx3      = tem
        tx4      = tx4 + eta(l) * hol(l)
        tx5      = tx5 + gms(l) * hol(l)
!
        qil(l)   = max(0.0, min(1.0, (tcr-tol(l))*tcrf))
        qll(l+1) = (0.5*alhf) * st2 * (qil(l)+qil(l+1)) + one
      enddo
!
!     for the cloud top -- l=kd
!
      l = kd
!
      tem      = qst(l) - gaf(l) * hst(l)
      tem1     = (tx3 + tem) * 0.5
      st2      = (gaf(l)+gaf(l+1)) * 0.5
!
      fco(l+1) =            tem1 + st2 * hbl
      rnn(l+1) = rnn(l+1) * tem1 + st2 * tx4
      gmh(l+1) = gmh(l+1) * tem1 + st2 * tx5
!
      fco(l)   = tem + gaf(l) * hbl
      rnn(l)   = tem * rnn(l) + (tx4 + eta(l)*hol(l)) * gaf(l)
      gmh(l)   = tem * gmh(l) + (tx5 + gms(l)*hol(l)) * gaf(l)
!
!   replace fco for the bottom
!
      fco(kbl) = qbl
      qis      = qii(kd)
      qls      = qli(kd)
      qil(kd)  =  max(0.0, min(1.0, (tcr-tol(kd))*tcrf))
      qll(kd1) = (0.5*alhf) * st2 * (qil(kd) + qil(kd1)) + one
      qll(kd ) = alhf * gaf(kd) * qil(kd) + one
!


      do l=kd,kb1
         fco(l) = fco(l+1) - fco(l)
         rnn(l) = rnn(l+1) - rnn(l) + eta(l)*qol(l)
         gmh(l) = gmh(l+1) - gmh(l) + gms(l)*qol(l)
!
!        tem    = c0 * eta(l)
         tem    = (c0*(1.0-qil(l)) + c0i*qil(l)) * eta(l)
         bkc(l) = qll(l+1) - tem * (1.0-akt(l))
         akt(l) = qll(l)   + tem * akt(l)
         akc(l) = 1.0 / akt(l)
      enddo

      rnn(kbl) = 0.0
      tx3      = 0.0
      tx4      = 0.0
      tx5      = 0.0
      do l=kb1,kd1,-1
        tem    = bkc(l-1)       * akc(l)
        tx3    = (tx3 + fco(l)) * tem
        tx4    = (tx4 + rnn(l)) * tem
        tx5    = (tx5 + gmh(l)) * tem
      enddo
      if (kd .lt. kb1) then
         st1   = tx3 / bkc(kd)
         hsd   = hst(kd1)
     *         + ltl(kd1) * (esfac *(qol(kd1)-qst(kd1))
     *                  +  st1 * (1.0 + esfac*qol(kd1)) )
     *         - (0.5*alhf) * st1 * (qil(kd1)+qil(kd))
      else
         hsd   = hbl
      endif
!
      tx3 = (tx3 + fco(kd)) * akc(kd)
      tx4 = (tx4 + rnn(kd)) * akc(kd)
      tx5 = (tx5 + gmh(kd)) * akc(kd)
      alm = alhf*qil(kd) - ltl(kd) * (1.0 + esfac*qol(kd))
!
      hsu = hst(kd) + ltl(kd) * esfac * (qol(kd)-qst(kd))
!
!===> vertical integrals needed to compute the entrainment parameter
!
      tx1 = alm * tx4
      tx2 = alm * tx5

      do l=kd,kb1
        tau = hol(l) - hsu
        tx1 = tx1 + tau * eta(l)
        tx2 = tx2 + tau * gms(l)
      enddo
!
!     modify hsu to include cloud liquid water and ice terms
!
      hsu   = hsu - alm * tx3
!
      clp   = zero
      alm   = -100.0
      hos   = hol(kd)
      qos   = qol(kd)
      unsat = hbl .gt. hsu .and. abs(tx1) .gt. 1.0e-4


!***********************************************************************


       st1    = half*(hsu + hsd)
       if (unsat) then
!
!  standard case:
!   cloud can be neutrally bouyant at middle of level kd w/ +ve lambda.
!   epp < .25 is required to have real roots.
!
         st1 = one / tx1
         x00 = -(hbl-hsu) * st1
         epp = -x00 * tx2 * st1
         akm = epp
         epp = epp + epp
!        if((tx1.lt.zero .and. epp.lt..5) .or.
!    *      (tx1.gt.zero .and. epp.lt.0.0)      ) then
!          if (epp.eq.zero) epp = 1.e-8
!          sgn = sign(1.,tx1)
!          alm = x00*(one/epp)*(one+sgn*sqrt(one-epp-epp))
!          clp = one
!        endif
!
         if (epp .lt. 0.5) then
           if (epp.eq.zero) epp = 1.e-8
           tem1   = x00 / epp
           st2    = tem1 * sqrt(one - epp - epp)
           tem2   = tem1 - st2
           tem1   = tem1 + st2
           if (tem1 .gt. almax) tem1 = -100.0
           if (tem2 .gt. almax) tem2 = -100.0
           alm    = max(tem1, tem2)
           clp    = one
         endif
!
!
!  clip case:
!   non-entrainig cloud detrains in lower half of top layer.
!   no clouds are allowed to detrain below the top layer.
!
       elseif ( (hbl .le. hsu) .and.
     *          (hbl .gt. st1   )     ) then
         alm = zero
         clp = (hbl-st1) / (hsu-st1)
       endif


!
      unsat = .true.
      if (almin1 .gt. 0.0) then
        if (alm .ge. almin1) unsat = .false.
      else
        lowest   = kd .eq. kb1
        if ( (alm .gt. zero) .or.
     *      (.not. lowest .and. alm .eq. zero) ) unsat = .false.
      endif
!
!===>  if no sounding meets second condition, return
!
      if (unsat) return 
!

      if(clp.gt.zero .and. clp.lt.one) then
        st1     = half*(one+clp)
        st2     = one - st1
        hst(kd) = hst(kd)*st1 + hst(kd1)*st2
        hos     = hol(kd)*st1 + hol(kd1)*st2
        qst(kd) = qst(kd)*st1 + qst(kd1)*st2
        qos     = qol(kd)*st1 + qol(kd1)*st2

        dlb(kd) = dlb(kd)*st1 + dlb(kd1)*st2
        ltl(kd) = ltl(kd)*st1 + ltl(kd1)*st2
        eta(kd) = eta(kd)*clp
        gms(kd) = gms(kd)*clp
      endif
!
      if (almin2 .ne. 0.0) then
        st1 = 0.0
        if (almin1 .ne. almin2) st1 = 1.0 / max(1.0e-10,(almin2-almin1))
          if (alm .lt. almin2) then
             clp = clp * (alm - almin1) * st1
          endif
      endif
!


!
!***********************************************************************
!
!    critical workfunction is included in this version
!
      if (crtfun) then
        acr = 0.0
        tem = prl(kd1) - (prl(kd1)-prl(kd)) * clp * half
        call crtwrk(tem, st1)
        acr = (prl(k) - tem) * st1
        clp = clp * qtl
      else
        acr = 0.0
        clp = clp * qtl
      endif
!
!===>  normalized massflux
!
!  eta is the thickness coming in and the mass flux going out.
!  gms is the thickness of the square; it is later reused for gamma_s
!
!     eta(k) = one

      do l=kb1,kd,-1
        eta(l) = eta(l+1) + alm * (eta(l) + alm * gms(l))
      enddo


!
!===>  cloud workfunction
!

      wfn   = zero
      akm   = zero
      det   = zero
      hcc   = hbl
      unsat = .false.
      qtl   = qst(kb1) - gaf(kb1)*hst(kb1)
      qtv   = qbl
      tx1   = hbl
!
      tem   = qst(kbl) - gaf(kbl)*hst(kbl)
      qtv   = 0.5 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
      det   = max(0.0, qbl-qtv)
      qtv   = qbl - det
!
      do l=kb1,kd1,-1
         del_eta = eta(l) - eta(l+1)
         hccp = hcc + del_eta*hol(l)
!

 
         qtlp = qst(l-1) - gaf(l-1)*hst(l-1)
         qtvp = 0.5 * ((qtlp+qtl)*eta(l)
     *              + (gaf(l)+gaf(l-1))*hccp)
         detp = (bkc(l)*det - (qtvp-qtv) + del_eta*qol(l)) * akc(l)
!
         tem1   = akt(l)   - qll(l)
         tem2   = qll(l+1) - bkc(l)
         rnn(l) = tem1*detp  + tem2*det

         if (detp .le. zero) unsat = .true.

         st1  = hst(l) - ltl(l)*esfac*(qst(l)-qol(l))


         tem2 = hccp   + detp   * (qil(l-1)+qil(l)) * (alhf*0.5)
!

         st2  = ltl(l) * (1.0 + esfac*qol(l))
         tem5 = qli(l) + qii(l)
         tem3 = (tx1  - eta(l+1)*st1 - st2*(det-tem5))  * dlb(l)
         tem4 = (tem2 - eta(l  )*st1 - st2*(detp-tem5)) * dlt(l)

         st1  = tem3 + tem4

         wfn = wfn + st1       
         akm = akm - min(st1,zero)
!

!        buy(l) = eta(l+1)*tem3 + eta(l)*tem4
         buy(l) = 0.5 * (eta(l+1) + eta(l)) * st1
!
         hcc = hccp
         det = detp
         qtl = qtlp
         qtv = qtvp
         tx1 = tem2


      enddo


      del_eta = eta(kd) - eta(kd1)
      hccp    = hcc + del_eta*hos
!

 
      qtlp = qst(kd) - gaf(kd)*hst(kd)
      qtvp = qtlp*eta(kd) + gaf(kd)*hccp
      detp = (bkc(kd)*det - (qtvp-qtv) + del_eta*qos) * akc(kd)
!
      tem1    = akt(kd)  - qll(kd)
      tem2    = qll(kd1) - bkc(kd)
      rnn(kd) = tem1*detp  + tem2*det
!
      if (detp.le.zero) unsat = .true.


!
!#ifdef DBG
      st1 = half * (hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
     *    + hst(kd1) - ltl(kd1)*esfac*(qst(kd1)-qol(kd1)))
      st2 = half * (ltl(kd)  * (1.0+esfac*qol(kd))
     *           +  ltl(kd1) * (1.0+esfac*qol(kd1)) )
!#endif
!     st1 = hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
!     st2 = ltl(kd)  * (1.0+esfac*qol(kd))
!
 

!     st1 = half*(tx1-eta(kd1)*st1-st2*det)*dlb(kd)
      tem5 = half * (qli(kd1)+qii(kd1)+qls+qis)
      st1  = half * (tx1-eta(kd1)*st1-st2*(det-tem5))*dlb(kd)
!

      wfn = wfn + st1
      akm = akm - min(st1,zero)
!

!     buy(kd) = eta(kd1) * st1
      buy(kd) = 0.5 * (eta(kd1) + eta(kd)) * st1
!
      det = detp
      hcc = hccp
      akm = akm / wfn


!***********************************************************************
!
!     if only to calculate workfunction save it and return
!
      if (wrkfun) then
        if (wfn .ge. 0.0) wfnc = wfn
        return
      elseif (.not. crtfun) then
        acr = wfnc
      endif
!
!===>  third check based on cloud workfunction
!
      calcup = .false.

      tem  =  min(cd*50.0, max_neg_bouy)
      if (wfn .gt. acr .and.  (.not. unsat)
     *                 .and. akm .le. tem) then
!


        calcup = .true.
      endif




!
!===>  if no sounding meets third condition, return
!

      if (.not. calcup) return
!


! 
      do l=kbl,k 
        rnn(l) = 0.0 
      enddo
!

!
!     if downdraft is to be invoked, do preliminary check to see
!     if enough rain is available and then call ddrft.
!
      ddft = .false.
      if (dndrft) then
!
        train = 0.0
        if (clp .gt. 0.0) then
          do l=kd,kb1
            train = train + rnn(l)
          enddo
        endif

        pl = (prl(kd1) + prl(kd))*half
        if (train .gt. 1.0e-4 .and. pl .le. pdd) ddft  = .true.
        do l=kd,km1
          if (l .lt. kbl-2 .and. ddft) then
            if (buy(l) .lt. 0.1) ddft = .false.
          endif
        enddo
      endif
!
      if (ddft) then
!
!     call downdraft scheme based on (cheng and arakawa, 1997)
!
        call ddrft(
     i              k, kd
     i,             cp, alhl, alhf, rkap, grav, tla, alfint
     i,             tol, qol, hol, prl, qst, hst, gam, gaf, hbl, qbl
     i,             qrb, qrt, buy, kbl, idh, eta, rnn
     *,             alm, wfn, train, ddft
     o,             etd, hod, qod, evp, dof, cldfrd, etz
     *,             gms, gsd, ghd)
!    *,             tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9)

      endif
!
!  no downdraft case (including case with no downdraft soln)
!
      if (.not. ddft) then
        do l=kd,k+1
          etd(l) = 0.0
          hod(l) = 0.0
          qod(l) = 0.0
        enddo
        do l=kd,k
          evp(l) = 0.0
          etz(l) = 0.0

        enddo

      endif
!

!
!===> calculate gammas  i.e. tendencies per unit cloud base massflux
!           includes downdraft terms!



      tx1 = qib
      tx2 = qlb
      do l=kb1,kd,-1
         del_eta = eta(l) - eta(l+1)
         tx1     = tx1 + del_eta * qii(l)
         tx2     = tx2 + del_eta * qli(l)
      enddo
!

!
      tem = det * qil(kd)
      st1 = (hcc+alhf*tem-eta(kd)*hst(kd)) / (1.0+gam(kd))
      ds  = eta(kd1) * (hos- hol(kd)
     *    - alhl*(qos - qol(kd)) + alhf*(qis-qii(kd)) )
      dh  = eta(kd1) * (hos- hol(kd))
      gms(kd) = (ds + st1) * pri(kd)
      gmh(kd) = pri(kd) * (hcc-eta(kd)*hos + dh)
!

!      tendency for suspended environmental ice and/or liquid water
!
      qil(kd) =     (tem + eta(kd1)*(qis-qii(kd))
     *                   + tx1 - eta(kd)*qis ) * pri(kd)
      qll(kd) = (det-tem + eta(kd1)*(qls-qli(kd))
     *                   + tx2 - eta(kd)*qls ) * pri(kd)
!
      ghd(kd) = 0.0
      gsd(kd) = 0.0
!
      do l=kd1,k
       st1 = one - alfint(l)
       st2 = one - alfinq(l)
       if (l .lt. kbl) then
         hb       = alfint(l)*hol(l-1) + st1*hol(l)
         qb       = alfint(l)*qol(l-1) + st1*qol(l)

         tem      = alfinq(l)*qii(l-1) + st2*qii(l)
         tem2     = alfinq(l)*qli(l-1) + st2*qli(l)
 
         tem1     = eta(l) * (tem - qii(l))
         tem3     = eta(l) * (tem2 - qli(l))

         tem5     =  etd(l) * (hod(l) - hb)
         tem6     =  etd(l) * (qod(l) - qb)
!
         dh       = eta(l) * (hb - hol(l)) + tem5
         ds       = dh - alhl * (eta(l) * (qb - qol(l)) + tem6)

         gmh(l)   = dh * pri(l)
         gms(l)   = ds * pri(l)
!
         ghd(l)   = tem5 * pri(l)
         gsd(l)   = (tem5 - alhl * tem6) * pri(l)
!


         qil(l)   = tem1 * pri(l)
         qll(l)   = tem3 * pri(l)

         tem1     = eta(l) * (qii(l-1) - tem)
         tem3     = eta(l) * (qli(l-1) - tem2)

         dh       = eta(l) * (hol(l-1) - hb) - tem5
         ds       = dh - alhl * eta(l) * (qol(l-1) - qb)
     *                 + alhl * (tem6 - evp(l-1))

         gmh(l-1) = gmh(l-1) + dh * pri(l-1)
         gms(l-1) = gms(l-1) + ds * pri(l-1)
!
         ghd(l-1) = ghd(l-1) - tem5 * pri(l-1)
         gsd(l-1) = gsd(l-1) - (tem5-alhl*(tem6-evp(l-1))) * pri(l-1)



         qil(l-1) = qil(l-1) + tem1 * pri(l-1)
         qll(l-1) = qll(l-1) + tem3 * pri(l-1)
       elseif (l .eq. kbl) then
         hb       = alfint(l)*hol(l-1) + st1*hbl
         qb       = alfint(l)*qol(l-1) + st1*qbl

!        hb       = hbl
!        qb       = qbl
         hbd      = alfint(l)*hol(l-1) + st1*hol(l)
         qbd      = alfint(l)*qol(l-1) + st1*qol(l)

         tem      = alfinq(l)*qii(l-1) + st2*qib
         tem2     = alfinq(l)*qli(l-1) + st2*qlb

         tem1     = eta(l) * (tem - qib)
         tem3     = eta(l) * (tem2 - qlb)

         tem5     =  etd(l) * (hod(l) - hbd)
         tem6     =  etd(l) * (qod(l) - qbd)

         tem4     = (grav/cmb2pa) * pris
         tx1      = eta(l) * (hb - hbl) * tem4
         tx2      = tx1 - alhl * eta(l) * (qb - qbl) * tem4
         dh       = tem5

         ds       =  dh - alhl * (tem6 + evp(l))


         gmh(l)   = tx1 + dh * pri(l)
         gms(l)   = tx2 + ds * pri(l)
!
         ghd(l)   = tem5 * pri(l)
         gsd(l)   = (tem5 - alhl * (tem6+evp(l))) * pri(l)
!


         qil(l)   = tem1 * tem4
         qll(l)   = tem3 * tem4

         tem1     = eta(l) * (qii(l-1) - tem)
         tem3     = eta(l) * (qli(l-1) - tem2)

         dh       = eta(l) * (hol(l-1) - hb) - tem5
         ds       = dh - alhl * eta(l) * (qol(l-1) - qb)
     *                 + alhl * (tem6 - evp(l-1))

         gmh(l-1) = gmh(l-1) + dh * pri(l-1)
         gms(l-1) = gms(l-1) + ds * pri(l-1)
!
         ghd(l-1) = ghd(l-1) - tem5 * pri(l-1)
         gsd(l-1) = gsd(l-1) - (tem5-alhl*(tem6-evp(l-1)))
     *                                  * pri(l-1)


         qil(l-1) = qil(l-1) + tem1 * pri(l-1)
         qll(l-1) = qll(l-1) + tem3 * pri(l-1)
       else
         hbd      = alfint(l)*hol(l-1) + st1*hol(l)
         qbd      = alfint(l)*qol(l-1) + st1*qol(l)
         tem5     =  etd(l) * (hod(l) - hbd)
         tem6     =  etd(l) * (qod(l) - qbd)
         dh       =  tem5
         ds       =  dh - alhl * (tem6 + evp(l))
!
         gmh(l)   = tx1 + dh * pri(l)
         gms(l)   = tx2 + ds * pri(l)
         ghd(l)   = dh * pri(l)
         gsd(l)   = ds * pri(l)
!
         dh       = - tem5
         ds       = dh  + alhl * tem6
         gmh(l-1) = gmh(l-1) + dh * pri(l-1)
         gms(l-1) = gms(l-1) + ds * pri(l-1)
!
         ghd(l-1) = ghd(l-1) + dh * pri(l-1)
         gsd(l-1) = gsd(l-1) + ds * pri(l-1)
!
         qil(l)   = qil(l-1)
         qll(l)   = qll(l-1)
       endif



      enddo

      hbd  = hol(k)
      qbd  = qol(k)
      tem5 =  etd(k) * (hod(k+1) - hbd)
      tem6 =  etd(k) * (qod(k+1) - qbd)
      dh   = - tem5
      ds   = dh  + alhl * tem6
      tem1 = dh * pri(k)
      tem2 = ds * pri(k)
      gmh(k) = gmh(k) + tem1
      gms(k) = gms(k) + tem2
      ghd(k) = ghd(k) + tem1
      gsd(k) = gsd(k) + tem2
!
      tem4   = - (grav/cmb2pa) * pris
      tx1    = dh * tem4
      tx2    = ds * tem4
!
      do l=kbl,k
        gmh(l) = gmh(l) + tx1
        gms(l) = gms(l) + tx2
        ghd(l) = ghd(l) + tx1
        gsd(l) = gsd(l) + tx2
      enddo
!


!***********************************************************************
!***********************************************************************

!===>  kernel (akm) calculation begins

!===>  modify sounding with unit mass flux
!

      testmb = 0.1

      do l=kd,k

         tem1   = gmh(l) - ghd(l)
         tem2   = gms(l) - gsd(l)
         hol(l) = hol(l) +  tem1*testmb
!        qol(l) = qol(l) + (tem1-tem2) * (testmb/alhl)
         qol(l) = qol(l) + (tem1-tem2
     *                       + alhf*qil(l))  * (testmb/alhl)
         hst(l) = hst(l) +  tem2*(one+gam(l))*testmb
         qst(l) = qst(l) +  tem2*gam(l)*(testmb/alhl)
      enddo
!

      hos = hos + (gmh(kd)-ghd(kd))  * testmb
      qos = qos + (gmh(kd)-gms(kd)
!    *          -  ghd(kd)+gsd(kd) ) * (testmb/alhl)
     *          -  ghd(kd)+gsd(kd)
     *          + alhf*qil(kd)) * (testmb/alhl)

!
      tem = prl(k+1) - prl(k)
      hbl = hol(k) * tem
      qbl = qol(k) * tem
      do l=km1,kbl,-1
        tem = prl(l+1) - prl(l)
        hbl = hbl + hol(l) * tem
        qbl = qbl + qol(l) * tem
      enddo
      hbl = hbl * pris
      qbl = qbl * pris


!***********************************************************************

!===>  cloud workfunction for modified sounding, then kernel (akm)
!
      akm = zero
      tx1 = zero
      qtl = qst(kb1) - gaf(kb1)*hst(kb1)
      qtv = qbl
      hcc = hbl
      tx2 = hcc
      tx4 = (alhf*0.5)*max(0.0,min(1.0,(tcr-tol(kb1))*tcrf))
!
      tem   = qst(kbl) - gaf(kbl)*hst(kbl)
      qtv   = 0.5 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
      tx1   = max(0.0, qbl-qtv)
      qtv   = qbl - tx1
!

      do l=kb1,kd1,-1
         del_eta = eta(l) - eta(l+1)
         hccp = hcc + del_eta*hol(l)
!
         qtlp = qst(l-1) - gaf(l-1)*hst(l-1)
         qtvp = 0.5 * ((qtlp+qtl)*eta(l)
     *                +(gaf(l)+gaf(l-1))*hccp)
         detp = (bkc(l)*tx1 - (qtvp-qtv) + del_eta*qol(l)) * akc(l)

         st1 = hst(l) - ltl(l)*esfac*(qst(l)-qol(l))

         tem2 = (alhf*0.5)*max(0.0,min(1.0,(tcr-tol(l-1))*tcrf))
         tem1 = hccp + detp * (tem2+tx4)



         st2  = ltl(l) * (1.0 + esfac * qol(l))
         tem5   = qli(l) + qii(l)
         akm  = akm +
     *     (  (tx2  -eta(l+1)*st1-st2*(tx1-tem5))  * dlb(l)
     *      + (tem1 -eta(l  )*st1-st2*(detp-tem5)) * dlt(l) )
!


         hcc  = hccp
         tx1  = detp
         tx2  = tem1
         qtl  = qtlp
         qtv  = qtvp
         tx4  = tem2
      enddo


!  eventhough we ignore the change in lambda, we still assume
!  that the cloud-top contribution is zero; as though we still
!  had non-bouyancy there.



!
!#ifdef DBG
      st1 = half * (hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
     *    +  hst(kd1) - ltl(kd1)*esfac*(qst(kd1)-qol(kd1))) 
      st2 = half * (ltl(kd)  * (1.0+esfac*qol(kd))
     *           +  ltl(kd1) * (1.0+esfac*qol(kd1)) )
!#endif
!
!     st1 = hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
!     st2 = ltl(kd)  * (1.0+esfac*qol(kd))
! 
      tem5 = half * (qli(kd1)+qii(kd1)+qls+qis)
      akm  = akm + half * (tx2-eta(kd1)*st1-st2*(tx1-tem5))*dlb(kd)

      akm = (akm - wfn) * (one/testmb)


!***********************************************************************

!===>   mass flux

      amb = - (wfn-acr) / akm

!===>   relaxation and clipping factors
!     tem    = max(min((prl(kd1)-prl(kd))*(rasalf*0.02), 0.5), 0.1)
!     amb(i) = amb * clp * tem
!

      amb = amb * clp * rasalf
       
!===>   sub-cloud layer depth limit on mass flux

      ambmax = (prl(kp1)-prl(kbl))*(cmb2pa*fracbl/grav)
      amb    = amax1(amin1(amb, ambmax),zero)



!***********************************************************************
!*************************results***************************************
!***********************************************************************





!===>  precipitation and clw detrainment


!
      dsfc = dsfc + amb * etd(k) * (1.0/dt)
!
      do l=kbl,kd,-1
          pcu(l) = pcu(l) + amb*rnn(l)      !  (a40)

      enddo


!
!===> temparature and q change and cloud mass flux due to cloud type kd
!
 

      tx1 = amb * (one/cp)
      tx2 = amb * (one/alhl)
      do l=kd,k
        st1    = gms(l)*tx1
        toi(l) = toi(l) + st1
        tcu(l) = tcu(l) + st1
        tcd(l) = tcd(l) + gsd(l) * tx1



        flx(l) = flx(l) + eta(l)*amb
!
        qii(l) = qii(l) + qil(l) * amb
        tem    = 0.0


        qli(l) = qli(l) + qll(l) * amb + tem



        st1          = (gmh(l)-gms(l)+alhf*qil(l))

     *                                   * tx2

        qoi(l) = qoi(l) + st1
        qcu(l) = qcu(l) + st1
        qcd(l) = qcd(l) + (ghd(l)-gsd(l)) * tx2



      enddo


!

      tx1 = 0.0
      tx2 = 0.0


!
!     reevaporation of falling convective rain
!
      if (revap) then
       gravcon = cmb2pa/grav
       afc     = -(1.04e-4*dt)*(3600./dt)**0.578
       rknob   = 5.0
!
       tx3    = amb*eta(kd)*pri(kd)

       cldfrd = min(amb*cldfrd, 1.0)
!

       do l=kd,k

!                                                 for l=kd,k
         if (l .ge. idh .and. ddft) then
           tx2 = tx2 + amb * rnn(l)
         else
           tx1 = tx1 + amb * rnn(l)
         endif
!

         clfrac = min(tx3*rknob, 1.0)


         if (tx1 .gt. 0. .or. tx2 .gt. 0.0) then
          teq     = toi(l)
          qeq     = qoi(l)
          pl      = 0.5 * (prl(l+1)+prl(l))

          st1     = max(0.0, min(1.0, (tcr-teq)*tcrf))
          st2     = st1*elfocp + (1.0-st1)*elocp

          call qsatcn ( teq,pl,qsteq,dqdt )
          deltaq = 0.5 * (qsteq-qeq) / (1.+st2*dqdt)
!
          qeq    = qeq + deltaq
          teq    = teq - deltaq*st2
!
          tem1   = max(0.0, min(1.0, (tcr-teq)*tcrf))
          tem2   = tem1*elfocp + (1.0-tem1)*elocp

          call qsatcn ( teq,pl,qsteq,dqdt )
          deltaq = (qsteq-qeq) / (1.+tem2*dqdt)
!
          qeq    = qeq + deltaq
          teq    = teq - deltaq*tem2

          if (qeq .gt. qoi(l)) then
            potevap = (qeq-qoi(l))*(prl(l+1)-prl(l))*gravcon

            tem3    = sqrt(pl*0.001)
            tem4    = potevap * (1. - exp( afc*sqrt(tx1*tem3) ) )
            actevap = min(tx1, tem4*clfrac)
!
            tem4    = potevap * (1. - exp( afc*sqrt(tx2*tem3) ) )
            tem4    = min(tx2, tem4*cldfrd)
!
            tx1     = tx1 - actevap
            tx2     = tx2 - tem4
            st1      = (actevap+tem4) * pri(l)
            qoi(l)  = qoi(l) + st1
            qcu(l)  = qcu(l) + st1
!

            st1     = st1 * elocp
            toi(l)  = toi(l) - st1 
            tcu(l)  = tcu(l) - st1
          endif
         endif
       enddo
!
      else
       do l=kd,k
         tx1 = tx1 + amb * rnn(l)

       enddo
      endif

      cup = cup + tx1 + tx2 + dof * amb


!
!    mixing of passive tracers
!
      do n=1,m

       do l=kd,k
         hol(l) = roi(l,n)
       enddo
!
        hcc     = rbl(n)
        hod(kd) = hol(kd)
!      compute downdraft properties for the tracer
       do l=kd1,k
         st1 = one - alfint(l)
         hb  = alfint(l)  * hol(l-1) + st1 * hol(l)
         if (etz(l-1) .ne. 0.0) then
           del_eta = etd(l) - etd(l-1)
           tem     = 1.0 / etz(l-1)
           if (del_eta .gt. 0.0) then
             hod(l) = (etd(l-1)*(hod(l-1)-hol(l-1))
     *                +  etd(l)  *(hol(l-1)-hb)
     *                +  etz(l-1)*hb) * tem
           else
             hod(l) = (etd(l-1)*(hod(l-1)-hb) + etz(l-1)*hb) * tem
           endif
         else
           hod(l) = hb
         endif
       enddo
             
       do l=kb1,kd,-1
          hcc = hcc + (eta(l)-eta(l+1))*hol(l)
       enddo
!
       gmh(kd) = pri(kd) * (hcc-eta(kd)*hol(kd))
       do l=kd1,k
        st1 = one - alfint(l)
         if (l .lt. kbl) then
           hb       = alfint(l) * hol(l-1) + st1 * hol(l)
           tem5     = etd(l)    * (hod(l) - hb)
           dh       = eta(l)    * (hb - hol(l)) + tem5
           gmh(l  ) = dh * pri(l)
           dh       = eta(l)    * (hol(l-1) - hb) - tem5
           gmh(l-1) = gmh(l-1)  + dh * pri(l-1)
         elseif (l .eq. kbl) then
           hb       = alfint(l) * hol(l-1) + st1 * rbl(n)
           hbd      = alfint(l) * hol(l-1) + st1 * hol(l)
           dh       = etd(l)    * (hod(l) - hbd)
           tem4     = (grav/cmb2pa) * pris
           tx1      = eta(l)    * (hb - rbl(n)) * tem4
           gmh(l)   = tx1       + dh * pri(l)
           dh       = eta(l)    * (hol(l-1) - hb) - dh
           gmh(l-1) = gmh(l-1)  + dh * pri(l-1)
         else
           hbd      = alfint(l) * hol(l-1) + st1 * hol(l)
           dh       = etd(l)    * (hod(l) - hbd)
           gmh(l)   = tx1       + dh * pri(l)
           gmh(l-1) = gmh(l-1)  - dh * pri(l-1)
         endif
       enddo
!
       do l=kd,k
         st1      = gmh(l)*amb
         roi(l,n) = hol(l)    + st1
         rcu(l,n) = rcu(l,n) + st1
       enddo
      enddo                             ! tracer loop m

!***********************************************************************
!***********************************************************************
!***********************************************************************

      return
      end
