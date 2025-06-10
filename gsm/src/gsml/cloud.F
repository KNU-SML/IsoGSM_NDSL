#include <define.h>
      subroutine cloud (len, lenc, k, km, dt, rasal, frac
     *,                 rhfacl, rhfacs
     *,                 cp,  alhl, grav, crmsfx, prj
     *,                 cd
c
     *,                 sgo, onbdsg, prh,  fpk, hpk, sge
     *,                 eps, del,    sdu,  sdd
     *,                 pii, pik, pki, gam
     *,                 uoi, voi, poi, qoi
     *,                 pcu, cln, tcu, qcu
#ifdef ISOTOPE
     *,                 qoi1,qoi2,pcui1,pcui2,ccndi1,ccndi2,ccnd
     *,                 eeqfo,eeqfd
#endif
     *,                 cwmdet,cmf,detr)
c   subroutine cloud
c        aug   2005   add detrained cloud water at the top of cumulus
c     16 march 1988
c
      parameter (zero=0.0, one=1.0, half=0.5, cmb2pa=100.0)
      parameter (ic=1,     daylen=86400.0,    epsln=0.001, qmin=1.0e-8)
      parameter (rhram=0.35)
      parameter(nx=7501)
#ifndef ICE
      dimension tbpvs0(nx)
      common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
c
c from mrf model
c
      parameter(rd= 2.8705e+2 , rv= 4.6150e+2
     &,         epsa=rd/rv,     epsam1=rd/rv-1.)
c
      parameter (zm1p04 = -1.04e-4, two89= 2.89e-5, zp44= 0.44)
      parameter ( zp01= 0.01, zp1 = 0.1, zp001= 0.001)
      parameter ( zp578 = 0.578, z3600 = 3600.)
      parameter ( thousand = 1000.)
c
      dimension gam(len,k)
      dimension pii(len),    pik(len),    pki(len)
      dimension uoi(len,k),  voi(len,k),  poi(len,km),  qoi(len,km)
     *,         cd(len)
#ifdef ISOTOPE
      dimension qoi1(len,km),qoi2(len,km),
     $          pcui1(len),pcui2(len),
     $          tmp(len,km),tmpi1(len,km),tmpi2(len,km),
     $          cond(len),condi1(len),condi2(len),
     $          cnd(len),cndi1(len),cndi2(len),
     $          ccnd(len),ccndi1(len),ccndi2(len)
      dimension eeqfo(len),eeqfd(len)
      integer   kk,flg
      real      ff,alpha1,alpha2,tmpi,alp1,alp2,evap,alp0,alp20
      real      freq !! function
      real      f1,f2,eqf
#ifdef NOFRAC
      parameter (f1=0.99D0,f2=0.60D0,eqf=0.0D0) !! f1:max rh, f2:min rh
#else
      parameter (f1=0.99D0,f2=0.60D0,eqf=0.45D0) !! f1:max rh, f2:min rh
#endif
      real      eqlee,eqlfr !! lee's equilibrium calculation
      parameter (eqlee=0,eqlfr=1.) !! 1:on, 0:off
      real      dif18o,difhdo   !! ratio of diffusion coef.
#ifdef NOFRAC
      parameter (dif18o=1.0D0,difhdo=1.0D0)
#else
      parameter (dif18o=1.02849D0,difhdo=1.02512D0)
#endif
      real      nn              !! degree of freedom
#ifdef NOFRAC
      parameter (nn=0.0)
#else
      parameter (nn=0.58)
#endif
      real      rh
      real      beta1,gamma1,rirat1,virat1,tmpi1o
      real      beta2,gamma2,rirat2,virat2,tmpi2o
      real      tmpq
#endif
c
      dimension tcu(len,km),  qcu(len,km),  cln(len,k), pcu(len)
c
      dimension onbdsg(km), fpk(km), hpk(km),  prh(km)
     *,         prj(km),    sgo(km), sge(km)
     *,         eps(km),    del(km), sdu(km),  sdd(km),  bta(km)
c
      real      cwmdet(len)
c
      real      cmf(len,km) 
      real      detr(len)    ! mass detrainment at cloud top. [-]
c
c     local
c
      real    hol(len ,k), qol(len,k), eta(len,k), hst(len,k),gmh(len,k)
      real    tx1(lenc), tx2(lenc), tx3(lenc), tx4(lenc), tx5(lenc)
      real    tx6(lenc), alm(lenc), wlq(lenc), wfn(lenc), akm(lenc)
      real    txc(lenc), txe(lenc), clf(lenc), uht(lenc), vht(lenc)
      real    qs1(lenc)
      integer ia(lenc),  i1(lenc),  i2(lenc)
      real    clp(lenc)
      real    cwmw(lenc)
#ifdef RASC2
      real    hu2qu(len,k),z(len,k)
#endif
c
c from qsatcn
c
      parameter ( h2omw  = 18.01, airmw  = 28.97)
      parameter ( esfac = h2omw/airmw       )
      parameter ( erfac = (1.0-esfac)/esfac )

      dimension tbsvp(nx), tbdsvp(nx)
      common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp
c
#ifdef CRAY_THREAD
cfpp$ expand (qsatcn, rncl, fpvs, fpvs0)
cfpp$ noconcur r
#endif
c
c
      km1 = k  - 1
      ic1 = ic + 1
      ip1 = 1
c
      onebg  = 1.0  / grav
      twobal = 2.0 / alhl
      albycp = alhl / cp
      el2orc = alhl*alhl / (rv*cp)
c
      do i=1,k
        bta(i) = fpk(i) * onebg
      enddo
      bta(ic)  = hpk(ic)  * onebg
c
      do l=ic,k
         do i=1,lenc
            tl = poi(i,l) * pik(i) * prh(l)
            pl = pii(i) * sgo(l)
c
            xj  = min(max(c1xsvp+c2xsvp*tl,1.),float(nx))
            jx  = min(xj,nx-1.)
            qx  = tbsvp(jx)  + (xj-jx)*(tbsvp(jx+1)-tbsvp(jx))
            dqx = tbdsvp(jx) + (xj-jx)*(tbdsvp(jx+1)-tbdsvp(jx))
            d = (pl-erfac*qx)
            if (d .gt. 0.) then
               d    = 1.0 / d
               gmh(i,l) = amin1(qx * d,1.0)
               dqs = (1.0 + erfac*gmh(i,l)) * d * dqx
            else
               gmh(i,l) = 1.0
               dqs = 0.
            endif
c
            gam(i,l) = 1.0 / (1.0 + dqs * albycp)
#ifdef RASC2
            hu2qu(i,l) =  dqs * albycp * gam(i,l) / alhl
#endif
         enddo
      enddo
c
      do i=1,lenc
        tem      = pik(i) * poi(i,k)
        qs1(i)   = gmh(i,k)
        qol(i,k) = amax1(amin1(qs1(i),qoi(i,k)), qmin)
        hol(i,k) = tem * (hpk(k)+cp*prh(k)) + qol(i,k)*alhl
        tx2(i)   = tem * fpk(k)
        eta(i,k) = zero
#ifdef RASC2
        z(i,k)   = zero
#endif
c
        clp(i)    = zero
        tem = rhfacs
c
        tem1   = qol(i,k) / qs1(i) - tem
        tx3(i) = tem1 + rhram
        if (tx3(i) .gt. 0.)
     *  clp(i) = max(zero, min(one, exp(20.0*tem1) )) ! 08/08/96
c
      enddo
c
      do l=km1,ic,-1
        do i=1,lenc
          tem      = pik(i) * poi(i,l)
          qs1(i)   = gmh(i,l)
          qol(i,l) = amax1(amin1(qs1(i),qoi(i,l)), qmin)
          tem1     = tx2(i) + tem * (hpk(l)+cp*prh(l))
          hol(i,l) = tem1   + qol(i,l) * alhl
          hst(i,l) = tem1   + qs1(i)   * alhl
          eta(i,l) = eta(i,l+1) + tem * bta(l)
          tx2(i)   = tx2(i)     + tem * fpk(l)
#ifdef RASC2
          z(i,l)   = eta(i,l)
#endif
        enddo
      enddo
c
      do i=1,lenc
        tx2(i) = hol(i,k)  - hst(i,ic)
        tx1(i) = zero
      enddo
c
c     entrainment parameter alm
c
      do l=ic,km1
        do i=1,lenc
          tx1(i) = tx1(i) + (hst(i,ic)-hol(i,l)) * (eta(i,l)-eta(i,l+1))
        enddo
      enddo
c
      len1 = 0
      len2 = 0
      isav = 0
      do i=1,lenc
        if (tx1(i) .gt. epsln .and. (tx2(i)  .gt. zero)
     *                        .and. (tx3(i)  .gt. zero) ) then
           len1      = len1 + 1
           ia(len1)  = i
           alm(len1) = tx2(i) / tx1(i)
        endif
      enddo
      len2 = len1
c
c     for non-entraining cloud lamda=0
c
      if (k .gt. 2) then
        do i=1,lenc
          if (tx2(i) .le. 0.0 .and. (hol(i,k) .gt. hst(i,ic1))
     *                        .and. (tx3(i)   .gt. zero) ) then
             len2      = len2 + 1
             ia(len2)  = i
             alm(len2) = 0.0
          endif
        enddo
      endif
c
      if (len2 .eq. 0) then
        return
      endif

      len11 = len1 + 1
c
c     normalized massflux
c
      do i=1,len1
        eta(i,k) = 1.0
        tx2(i) = pii(ia(i)) * sgo(ic)
        txc(i) = 1.0
      enddo
c
      do i=len11,len2
        eta(i,k) = 1.0
        wfn(i)   = 0.0
        ii       = ia(i)
        if (hst(ii,ic1) .lt. hst(ii,ic)) then
           txc(i) = (hst(ii,ic1)-hol(ii,k))/(hst(ii,ic1)-hst(ii,ic))
        else
           txc(i) = 0.0
        endif
        tx2(i) = pii(ii) * (sgo(ic1) - txc(i) * (sgo(ic1)-sgo(ic)))
      enddo
c
      call acritn(len2, tx2, tx3)
c
      do l=km1,ic,-1
        do i=1,len2
          tx5(i) = 1.0 + alm(i) * eta(ia(i),l)
        enddo
        do i=1,len2
          eta(i,l) = tx5(i)
        enddo
      enddo
c
         do i=1,len2
            tx5(i) = 1.0
         enddo
c
c     cloud workfunction
c
      if (len1 .gt. 0) then
        do i=1,len1
         wfn(i) = - gam(ia(i),ic) * eps(ic) * hst(ia(i),ic) * eta(i,ic1)
        enddo
      endif
c
      do i=1,len2
        tx1(i) = hol(ia(i),k)
        tem    = sgo(ic1) - txc(i) * (sgo(ic1)-sgo(ic))
        tx3(i) = tx3(i) * pii(ia(i)) * (sge(km1)-tem)
        akm(i) = 0.0
      enddo
c
      if (ic1 .le. km1) then
         do 380 l=km1,ic1,-1
         do 380 i=1,len2
         tem    = tx1(i) + (eta(i,l) - eta(i,l+1)) * hol(ia(i),l)
c
         wfnt   = gam(ia(i),l) * (tx1(i)*eps(l) + tem*del(l)
     *          - (eta(i,l+1)*eps(l) + eta(i,l)*del(l)) * hst(ia(i),l))
         wfn(i) = wfn(i) + wfnt
         akm(i) = akm(i) - min(wfnt, 0.0)
c
         tx1(i) = tem
  380    continue
      endif
c
      if (len1 .gt. 0) then
        do i=1,len1
          wfn(i) = wfn(i) + tx1(i) * gam(ia(i),ic) * eps(ic)
        enddo
      endif
      do i=1,len2
         if (wfn(i) .gt. 0.0) then
           akm(i) = akm(i) / wfn(i)
           tem    =  min(cd(ia(i))*100.00,0.50)
           if (akm(i) .gt. tem) wfn(i) = 0.0
         endif
         wfn(i) = wfn(i) - tx3(i)
      enddo

      lena = 0
      if (len1 .gt. 0) then
        do i=1,len1
          if (wfn(i) .gt. 0.0 .and. tx5(i) .ge. 0.0) then
            lena = lena + 1
            i1(lena) = ia(i)
            i2(lena) = i
            tx1(lena) = wfn(i)
            tx2(lena) = qs1(ia(i))
            txc(lena) = 1.0
          endif
        enddo
      endif
      lenb = lena
      do i=len11,len2
        if (wfn(i) .gt. 0.0 .and. txc(i) .gt. 0.0
     *                      .and. tx5(i) .ge. 0.0) then
           lenb = lenb + 1
           i1(lenb)  = ia(i)
           i2(lenb)  = i
           tx1(lenb) = wfn(i)
           tx2(lenb) = qs1(ia(i))
           tx4(lenb) = txc(i)
        endif
      enddo
c
      if (lenb .eq. 0) then
        return
      endif
c
      do i=1,lenb
        wfn(i) = tx1(i)
        qs1(i) = tx2(i)
      enddo
c
      do l=ic,k
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
        do i=1,lenb
          eta(i,l) = eta(i2(i),l)
        enddo
      enddo
c
      lena1 = lena + 1
c
      do i=1,lena
        ii = i1(i)
        tx6(i) = hst(ii,ic) - hol(ii,ic)
        uht(i) = uoi(ii,ic)
        vht(i) = voi(ii,ic)
      enddo
      do i=lena1,lenb
        ii = i1(i)
        txc(i) = tx4(i)
        tem    = txc(i) * (hol(ii,ic)-hol(ii,ic1)) + hol(ii,ic1)
        tx6(i) = hol(ii,k) - tem
c
        tem1   = txc(i) * (qol(ii,ic)-qol(ii,ic1)) + qol(ii,ic1)
        tx5(i) = tem    - tem1 * alhl
        qs1(i) = tem1   + tx6(i)*(one/alhl)
        tx3(i) = hol(ii,ic)
c
        uht(i) = txc(i) * (uoi(ii,ic)-uoi(ii,ic1)) + uoi(ii,ic1)
        vht(i) = txc(i) * (voi(ii,ic)-voi(ii,ic1)) + voi(ii,ic1)
      enddo
c
      c0 = 0.002              ! auto conversion coefficient

      do i=1,lenb
        ii = i1(i)
        wlq(i) = qol(ii,k)

#ifdef RASC2
        tx1(i) = hol(ii,k)   ! h*eta        
        tx4(i) = zero         ! sum of rain
#endif

        uht(i) = uoi(ii,k) - uht(i) * eta(i,ic)
        vht(i) = voi(ii,k) - vht(i) * eta(i,ic)
        txe(i) = hol(ii,k)
      enddo
c
      do l=km1,ic,-1
        do i=1,lenb
          ii = i1(i)
          tem    = eta(i,l) - eta(i,l+1)    ! d(eta)
          wlq(i) = wlq(i) + tem * qol(ii,l) ! (qu+qc)*eta

#ifdef RASC2
          tx1(i) = tx1(i) + tem * hol(ii,l) ! h*eta
          hu = tx1(i) / eta(i,l)                 ! hu
          qu = max( 0.0, gmh(ii,l) + hu2qu(ii,l) * (hu - hst(ii,l)) ) ! qu
          qceta = max( 0.0, wlq(i) - qu*eta(i,l) )                    ! qc*eta
          rainil = c0 * qceta * (z(ii,l) - z(ii,l+1))   ! rain at this level
          rainil = min( rainil, qceta )
          wlq(i) = wlq(i) - rainil        ! (qu+qc)*eta
          tx4(i) = tx4(i) + rainil        ! sum of rain
#endif

          uht(i) = uht(i) + tem * uoi(ii,l)
          vht(i) = vht(i) + tem * voi(ii,l)
        enddo
      enddo

      do i=1,lenb
        ii = i1(i)
        wlq(i) = wlq(i) - qs1(i) * eta(i,ic)  ! qc*eta
c---------------------------------------------------
c       when wlq<0 -> all cloud water becomes rain
c---------------------------------------------------
        if (wlq(i).lt.0.0) then
          tx4(i) = tx4(i) + wlq(i)  ! sum of rain
          wlq(i) = 0.0              ! detarained ... = 0
        endif
      end do
c
c     calculate gs and part of akm (that requires eta)
c
      do i=1,lenb
        ii = i1(i)
        tem        = (poi(ii,km1) - poi(ii,k)) * pik(ii)
        hol(i,k)   = tem    * (sdu(k)*onbdsg(k))
        hol(i,km1) = tem    * (sdd(km1)*onbdsg(km1))

        akm(i)     = zero
        tx2(i)     = pii(ii)*sgo(ic)
      enddo
c

      if (ic1 .le. km1) then
         do 750 l=km1,ic1,-1
         do 750 i=1,lenb
         ii = i1(i)
         tem        = (poi(ii,l-1) - poi(ii,l)) * pik(ii) * eta(i,l)
c
         hol(i,l)   = tem * (sdu(l  )*onbdsg(l))  + hol(i,l)
         hol(i,l-1) = tem * (sdd(l-1)*onbdsg(l-1))
c
         akm(i) = akm(i) - hol(i,l)*(eta(i,l)*del(l)+eta(i,l+1)*eps(l))
  750    continue
c
      endif
c
      call rncl(lenb, tx2, tx1, clf)
c
      do i=1,lenb

#ifdef RASC2
        tx2(i)  = zero                          ! evaporation
        cwmw(i) = wlq(i)                        ! detrained cloud water at cloud top(i=ic)
        wlq(i)  = tx4(i)                        ! precipitation 
#else

#ifdef RASC1
        tx2(i)  = zero                          ! evaporation
        cwmw(i) = (one - tx1(i)) * wlq(i)       ! detrained cloud water at cloud top(i=ic)
#else
        tx2(i) = (one - tx1(i)) * wlq(i)        ! evaporation
        cwmw(i) = zero                          ! detrained cloud water at cloud top(i=ic)
#endif
        wlq(i) = wlq(i) - tx2(i) - cwmw(i)      ! precipitation 

#endif

        tx1(i) = hol(i,ic)
      enddo
c
      do i=lena1,lenb
        ii = i1(i)
        tx1(i) = tx1(i) + (tx5(i)-tx3(i)+qol(ii,ic)*alhl) * onbdsg(ic)
      enddo

      do i=1,lenb
        hol(i,ic) = tx1(i) - tx2(i) * (alhl*onbdsg(ic))
      enddo
c
      if (lena .gt. 0) then
        do i=1,lena
          akm(i) = akm(i) - eta(i,ic1) * eps(ic) * tx1(i)
        enddo
      endif
c
      temu = (prh(k) - prj(k)) / (prh(k) - prh(k-1))
      temd = 1.0 - temu
c
      do i=1,lenb
        ii = i1(i)
        tx3(i)   =  temu*qol(ii,km1) + temd*qol(ii,k)
        gmh(i,k) = hol(i,k) + (tx3(i)-qol(ii,k)) * (alhl*onbdsg(k))
        akm(i) = akm(i) + gam(ii,km1)*eps(km1) * gmh(i,k)
      enddo
c
      if (ic1 .le. km1) then
         do 840 l=km1,ic1,-1
         temu = (prh(l) - prj(l)) / (prh(l) - prh(l-1))
         temd = 1.0 - temu
c
         do 840 i=1,lenb
         ii = i1(i)
         tx2(i) = tx3(i)
         tx3(i)   =  temu*qol(ii,l-1) + temd*qol(ii,l)
c
         gmh(i,l) = hol(i,l) + ((qol(ii,l) - tx2(i)) * eta(i,l+1)
     *                       +  (tx3(i)    - qol(ii,l)) * eta(i,l))
     *                       *  (alhl*onbdsg(l))
  840    continue
      endif
c
      do i=lena1,lenb
         tx3(i) = (txe(i)-tx6(i)-tx5(i)) - tx3(i) * alhl
      enddo
      do i=1,lena
         tx3(i) = (qol(i1(i),ic) - tx3(i)) * alhl
      enddo
      do i=1,lenb
         gmh(i,ic) = tx1(i) + onbdsg(ic)
     *                      * (eta(i,ic1)*tx3(i) + eta(i,ic)*tx6(i))
      enddo
c
c     calculate hc part of akm
c
      if (ic1 .le. km1) then
        do i=1,lenb
          tx1(i) = gmh(i,k)
        enddo
        do l=km1,ic1,-1
          do i=1,lenb
            ii = i1(i)
            tx1(i) = tx1(i) + (eta(i,l) - eta(i,l+1)) * gmh(i,l)
            tx2(i) = gam(ii,l-1) * eps(l-1)
          enddo
c
          if (l .eq. ic1) then
            do i=lena1,lenb
              tx2(i) = zero
            enddo
          endif
c
          do i=1,lenb
            akm(i) = akm(i) + (tx2(i) + gam(i1(i),l)*del(l)) * tx1(i)
          enddo
        enddo
      endif
c
      do i=lena1,lenb
        tx1(i) = sgo(ic) + (sgo(ic1) - sgo(ic)) * (one-txc(i))
        if ((tx1(i) .ge. sge(ic)) .and. (tx1(i) .le. sgo(ic1))) then
         txc(i)     = one - (tx1(i) - sge(ic)) / (sgo(ic1) - sge(ic))
c
         hol(i,ic1) = hol(i,ic1) + hol(i,ic) * (onbdsg(ic1)/onbdsg(ic))
         hol(i,ic)  = zero
c
         gmh(i,ic1) = gmh(i,ic1) + gmh(i,ic)*(onbdsg(ic1)/onbdsg(ic))
         gmh(i,ic)  = zero
        elseif (tx1(i) .lt. sge(ic)) then
         txc(i) = 1.0
        else
         txc(i) = 0.0
        endif
      enddo
c
      tem0 = pcu(i1(ip1))
      do i=1,lenb
        ii = i1(i)
        tem = wlq(i) + cwmw(i) ! precipitation + detrained cloud water
        if (akm(i) .lt. zero .and. tem .ge. 0.0) then
           wfn(i) = - clp(ii) * txc(i) * wfn(i) * rasal / akm(i)
        else
           wfn(i) = zero
        endif
        wfn(i) = amin1(wfn(i), frac)
c
c     compute cloud amount
c
        if (txc(i) .gt. 0.9999) then
           tx1(i) = cln(ii,ic)
        else
           tx1(i) = cln(ii,ic1)
        endif
        if (wfn(i) .gt. crmsfx)  tx1(i) = tx1(i) + clf(i)
        if (tx1(i) .gt. one)  tx1(i) = one
c
c     precipitation
c
        wlq(i) = (cmb2pa/grav) * wlq(i) * wfn(i) * pii(ii)  ! [mm]=[kg/m2]
c
c     detarained cloud water at cloud top
c
        cwmw(i) = cwmw(i) * wfn(i) * onbdsg(ic)    ! [kg/kg]
      enddo
c
      do i=1,lenb
        ii = i1(i)
        if (txc(i) .gt. 0.9999) then
           cln(ii,ic) = tx1(i)
        else
           cln(ii,ic1) = tx1(i)
        endif
      enddo
c
c     theta and q change due to cloud type ic
c
#ifdef ISOTOPE
      do i=1,lenb
        ii=i1(i)
        txe(i)=(grav/cmb2pa)/pii(ii)
        do l=k,ic,-1
          tmp(i,l)=0.D0
          tmpi1(i,l)=0.D0
          tmpi2(i,l)=0.D0
        enddo
      enddo
c
      do l=k,ic,-1
#else
      do l=ic,k
#endif
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
        do i=1,lenb
          ii = i1(i)
#ifdef ISOTOPE
          tmpq=qoi(ii,l)
#endif
c
          tem      = (gmh(i,l) - hol(i,l)) * wfn(i)
          tem1     =  hol(i,l) * wfn(i)
          poi(ii,l) = poi(ii,l) + tem1 * pki(ii) * (one/(cp*prh(l)))
          qoi(ii,l) = qoi(ii,l) + tem * (one/alhl)
          tcu(ii,l) = tcu(ii,l) + tem1 * (daylen/(cp*dt))
          qcu(ii,l) = qcu(ii,l) + tem  * (daylen/(alhl*dt))
          cmf(ii,l) = cmf(ii,l) + wfn(i)*eta(i,l) ! M=Mb*eta
#ifdef ISOTOPE
          tmp(i,l)=-tem*(one/alhl)/txe(i)/onbdsg(l)
          pl = pii(ii) * sgo(l) * cmb2pa
          tl = poi(ii,l) * pik(ii) * prh(l) 
          alpha1=freq(1,tl)
          alpha2=freq(2,tl)
#ifdef ICE
          es   = 10.0*fpvs(tl)
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tl,1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          es   = 10.0*fpvs01
#endif
          qs   = epsa * es / (tl + epsam1*es)
          rh=0.95  !! because above equation gives very low value for RH.
          if (tmpq.gt.0.D0) then
             virat1=qoi1(ii,l)/tmpq !! current vap isotop ratio
             virat2=qoi2(ii,l)/tmpq
          else
             virat1=0.
             virat2=0.
          endif
          if (tl.gt.273.15) then 
c for simple Rayleigh, force rh=0 and nn=0
            if ((virat1.gt.0.D0).and.(virat2.gt.0.D0)) then              
          gamma1=alpha1*rh/
     $         (alpha1*dif18o**nn*(rh-1)+1.D0)
          beta1=1.D0/(alpha1*dif18o**nn*(1-rh))-1.D0
          gamma2=alpha2*rh/
     $         (alpha2*difhdo**nn*(rh-1)+1)
          beta2=1.D0/(alpha2*difhdo**nn*(1-rh))-1.D0
            else
          gamma1=1.D0
          beta1=0.D0
          gamma2=1.D0
          beta2=0.D0
            endif              
          else !! ice -> vap
          gamma1=1.D0
          beta1=0.D0
          gamma2=1.D0
          beta2=0.D0
          endif
*** Super saturation Kinetic fractionation ***
* alpha_eff=alpha_eq*alpha_kin
* alpha_kin=S/(alpha_eq*D_iso/D*(S-1)+1)
* S=1-0.003*T [C] (Jouzel and Merlivat 1984)
          if (tl.le.253.15) then
            alpha1=alpha1*
     $            (1-0.003*(tl-273.15))/
     $            (alpha1*dif18o*(-0.003*(tl-273.15))+1)
            alpha2=alpha2*
     $            (1-0.003*(tl-273.15))/
     $            (alpha2*difhdo*(-0.003*(tl-273.15))+1)
*** linear interpolation between -20 to 0 C ***
          elseif (tl.le.273.15) then
            alp0=freq(1,273.15)
            alp20=freq(1,253.15)*
     $           (1-0.003*(-20.))/
     $           (freq(1,253.15)*dif18o*(-0.003*(-20.))+1)
            alpha1=
     $           (alp0*(tl-253.15)
     $           +alp20*(273.15-tl))/20.
            alp0=freq(2,273.15)
            alp20=freq(2,253.15)*
     $           (1-0.003*(-20.))/
     $           (freq(2,253.15)*difhdo*(-0.003*(-20.))+1)
            alpha2=
     $           (alp0*(tl-253.15)
     $           +alp20*(273.15-tl))/20.
          endif
          flg=0
c
          if (tem.lt.0.) then   !! condensation ( q > q' )
c
          if (qoi(ii,l).ge.0.D0) then !! q > q' >= 0.
          flg=1
c
          if (qoi(ii,l)-tem*(one/alhl).gt.0.D0) then
            ff=qoi(ii,l)/(qoi(ii,l)-tem*(one/alhl))
          else
            print*,'Something wrong1',i,l,qoi(ii,l)-tem*(one/alhl)
            ff=1.D0
          endif
          tmpi1(i,l)=tmpi1(i,l)+
     $         qoi1(ii,l)*(1.D0-ff**alpha1)
     $         /txe(i)/onbdsg(l)
          qoi1(ii,l)=qoi1(ii,l)*ff**alpha1
          tmpi2(i,l)=tmpi2(i,l)+
     $         qoi2(ii,l)*(1.D0-ff**alpha2)
     $         /txe(i)/onbdsg(l)
          qoi2(ii,l)=qoi2(ii,l)*ff**alpha2
c
          elseif (qoi(ii,l).lt.tem*(one/alhl)) then !! -q > -q' > 0.
          flg=2
c
          if ((qoi(ii,l)-tem*(one/alhl).ne.0.D0).and.
     $         (qoi(ii,l).ne.0.D0) ) then
            ff=(qoi(ii,l)-tem*(one/alhl))/qoi(ii,l)
          else
            print*,'Something wrong2',i,l,qoi(ii,l)-tem*(one/alhl)
            ff=1.D0
          endif
c
         alpha1=1.D0
         alpha2=1.D0
c          
          tmpi1(i,l)=tmpi1(i,l)+
     $         qoi1(ii,l)*(1.D0-ff**-alpha1)
     $         /txe(i)/onbdsg(l)
          qoi1(ii,l)=qoi1(ii,l)*ff**-alpha1
          tmpi2(i,l)=tmpi2(i,l)+
     $         qoi2(ii,l)*(1.D0-ff**-alpha2)
     $         /txe(i)/onbdsg(l)
          qoi2(ii,l)=qoi2(ii,l)*ff**-alpha2
c
          else !! q >= 0.  q' =< 0
          flg=3
c
          if (((tmp(i,l)+qoi(ii,l)/txe(i)/onbdsg(l)).gt.0.D0)
     $         .and.(tmp(i,l).ne.0.D0 ) )then
            ff=(tmp(i,l)+qoi(ii,l)/txe(i)/onbdsg(l))/tmp(i,l)
          else
            print*,'Something wrong3',i,l,tmp(i,l)
            ff=1.D0
          endif
c
          alpha1=1.D0
          alpha2=1.D0
c          
          tmpi=qoi1(ii,l)
          qoi1(ii,l)=
     $         (tmpi1(i,l)+tmpi/txe(i)/onbdsg(l))*
     $         (1.D0-ff**(-1./alpha1))*txe(i)*onbdsg(l)
          tmpi1(i,l)=(tmpi1(i,l)+tmpi/txe(i)/onbdsg(l))*ff**(-1./alpha1)
          tmpi=qoi2(ii,l)
          qoi2(ii,l)=
     $         (tmpi2(i,l)+tmpi/txe(i)/onbdsg(l))*
     $         (1.D0-ff**(-1./alpha2))*txe(i)*onbdsg(l)
          tmpi2(i,l)=(tmpi2(i,l)+tmpi/txe(i)/onbdsg(l))*ff**(-1./alpha2)
c
          endif
c
          elseif (tem.gt.0.D0) then !! evaporation
c     
          tmp(i,l)=0.D0
          evap=tem*(one/alhl)/txe(i)/onbdsg(l)
          do ll=l+1,k
            tmpi=tmp(i,ll)
            tmp(i,ll)=max(0.D0,tmp(i,ll)-evap)
            evap=max(0.D0,evap-(tmpi-tmp(i,ll)))
            if (tmpi.gt.0.D0) then
              ff=tmp(i,ll)/tmpi
c
c
ccc Currently, Stewart is very unstable for uplifting evaporation ccc
*** Stewart 1975 ***
*   rni/rn=(rni/ro-gamma*qi/q)*ff**beta+gamma*qi/q
*   gamma=alpha_eq*RH/(alpha_eq*(D_iso/D)**n*(RH-1)+1)
*   beta=1/(alpha_eq*(D_iso/D)**n*(1-RH))-1
*   RH: relative humidity
c
              rirat1=tmpi1(i,ll)/tmpi
              rirat2=tmpi2(i,ll)/tmpi
              tmpi1o=tmpi1(i,ll)
              tmpi2o=tmpi2(i,ll)
c              if ((rirat1.gt.0.D0).and.(rirat2.gt.0.D0)
c     $             .and.(virat1.gt.0.D0).and.(virat2.gt.0.D0)) then
c              tmpi1(i,ll)=
c     $             (rirat1*(eqf*ff**beta1+1.D0-eqf)
c     $             +eqf*(1.D0-ff**beta1)*gamma1*virat1)
c     $             *tmp(i,ll)
c              tmpi2(i,ll)=
c     $             (rirat2*(eqf*ff**beta2+1.D0-eqf)
c     $             +eqf*(1.D0-ff**beta2)*gamma2*virat2)
c     $             *tmp(i,ll)
c              else
              tmpi1(i,ll)=ff*tmpi1(i,ll)
              tmpi2(i,ll)=ff*tmpi2(i,ll)
c              endif              
              qoi1(ii,l)=
     $             qoi1(ii,l)+(tmpi1o-tmpi1(i,ll))*txe(i)*onbdsg(l)
              qoi2(ii,l)=
     $             qoi2(ii,l)+(tmpi2o-tmpi2(i,ll))*txe(i)*onbdsg(l)
c
c              if (tl.le.273.15) then !! ice -> vap
c                 alpha1=1.D0
c                 alpha2=1.D0
c              endif
c              alp1=(alpha1-1.D0)*eqf+1.D0
c              alp2=(alpha2-1.D0)*eqf+1.D0
c              if (ff.le.f1) then
c                 alp1=(alp1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c                 alp2=(alp2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c              endif
c
c              qoi1(ii,l)=qoi1(ii,l)+
c     $             tmpi1(i,ll)*(1.D0-ff**(1./alp1))*txe(i)*onbdsg(l)
c              tmpi1(i,ll)=tmpi1(i,ll)*ff**(1./alp1)
c              qoi2(ii,l)=qoi2(ii,l)+
c     $             tmpi2(i,ll)*(1.D0-ff**(1./alp2))*txe(i)*onbdsg(l)
c              tmpi2(i,ll)=tmpi2(i,ll)*ff**(1./alp2)
c
            else
              if ((tmpi.lt.0.D0).or.(tmpi1(i,ll).ne.0.D0)) then
                print*, i,l,ll,"strange1",tmpi,tmpi1(i,ll)
              endif
            endif
            if (evap.le.1.D-15) then
              flg=3+ll-l
              goto 100
            endif
          enddo
 100      continue
c
          endif
*** end for isotope updraft ***
c
#endif
        enddo
c
      enddo
      do i=1,lenb
        ii = i1(i)
        detr(ii) = detr(ii) + cmf(ii,ic) * onbdsg(ic)    ! mass detrainment at clodu top. [-]
      end do
      do i=1,lenb
        ii = i1(i)
c
c       cumulus precipitation
c
c       cumulus friction at the bottom level
c
        tem      = wfn(i) * (0.5*onbdsg(k))
        tx1(i)   = uoi(ii,km1) - uoi(ii,k)
        tx2(i)   = voi(ii,km1) - voi(ii,k)
c
        uoi(ii,k) = uoi(ii,k) + tx1(i) * tem
        voi(ii,k) = voi(ii,k) + tx2(i) * tem
      enddo
c
c     cumulus friction at all other levels
c
      do l=km1,ic1,-1
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
        do i=1,lenb
          ii = i1(i)
          tem = wfn(i) * (0.5*onbdsg(l))
          tem1   = tx1(i)
          tem2   = tx2(i)
          tx1(i) = (uoi(ii,l-1) - uoi(ii,l)) * eta(i,l)
          tx2(i) = (voi(ii,l-1) - voi(ii,l)) * eta(i,l)
c
          uoi(ii,l) = uoi(ii,l) + (tx1(i) + tem1) * tem
          voi(ii,l) = voi(ii,l) + (tx2(i) + tem2) * tem
        enddo
      enddo
c
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
      do i=1,lenb
        ii = i1(i)
        if (txc(i) .ge. 1.0) then
          tem    = wfn(i) * (onbdsg(ic)*0.5)
          uoi(ii,ic) = uoi(ii,ic) + (tx1(i) + uht(i) + uht(i)) * tem
          voi(ii,ic) = voi(ii,ic) + (tx2(i) + vht(i) + vht(i)) * tem
        elseif (txc(i) .gt. 0.0) then
          tem    = wfn(i) * (onbdsg(ic1)*0.5)
          uoi(ii,ic1) = uoi(ii,ic1) - (tx1(i) - uht(i) + uht(i)) * tem
          voi(ii,ic1) = voi(ii,ic1) - (tx2(i) - vht(i) + vht(i)) * tem
        endif
      enddo
c
      ipa = 1
      if (ipa .gt. lenb) ipa = 1
      ipo = i1(ipa)
c
      rphf = z3600/dt
      tem = onbdsg(ic)
      do i=1,lenb
        tx5(i) = tem * wfn(i) * eta(i,ic)
        if (tx5(i) .gt. 1.0) tx5(i) = 1.0
      enddo
c
       do  i = 1,lenb
         ii = i1(i)
         tx6(i) = wlq(i)
         tx4(i) = sqrt(pii(ii))
         txe(i) = (grav/cmb2pa) / pii(ii)
#ifdef ISOTOPE
*** isotope (temporary balancing) ***
         cond(i)=0.D0
         do l=ic,k
           cond(i)=cond(i)+tmp(i,l)
         enddo
         if (cond(i).ne.0.D0) then
           do l=ic,k
             tmpi1(i,l) = tmpi1(i,l)*wlq(i)/cond(i)
             tmpi2(i,l) = tmpi2(i,l)*wlq(i)/cond(i)
             tmp(i,l) = tmp(i,l)*wlq(i)/cond(i)
           enddo
         else
           if (wlq(i).gt.0.D0) then
             print*,'Something wrong',wlq(i),tmp(i,l)
           endif
         endif
         cond(i)=0.D0
         condi1(i)=0.D0
         condi2(i)=0.D0
c
         cnd(i)=0.
         cndi1(i)=0.
         cndi2(i)=0.
         do l=ic,k
           cnd(i)=cnd(i)+tmp(i,l)
           cndi1(i)=cndi1(i)+tmpi1(i,l)
           cndi2(i)=cndi2(i)+tmpi2(i,l)
         enddo
#endif
       enddo
c
c      do loop for moisture evaporation ability and convec evaporation.
c
       do 5000 l=ic,km
         tem    = rphf * sqrt(sgo(l)*zp001)
         dsgbg  = (onebg*cmb2pa) / onbdsg(l)
         rknobt = 5.0
c
         albyp = albycp / prh(l)
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
       do i=1,lenb
          ii = i1(i)
          tx1(i) = poi(ii,l) * pik(ii) * prh(l)
          tx3(i) = qoi(ii,l)
          txc(i) = sgo(l)*pii(ii)
c
c   first iteration - increment halved
c
#ifdef ICE
          es   = 10.0*fpvs(tx1(i))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          es   = 10.0*fpvs01
c         es   = 10.0*fpvs0(tx1(i))
#endif
          qs   = epsa * es / (txc(i) + epsam1*es)
          tsq  = tx1(i) * tx1(i)
          delq = 0.5 * (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
c
          tx2(i) = delq
          tx1(i) = tx1(i) - delq * albycp
          tx3(i) = tx3(i) + delq
c
c   second iteration
c
#ifdef ICE
          es   = 10.0*fpvs(tx1(i))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          es   = 10.0*fpvs01
c         es   = 10.0*fpvs0(tx1(i))
#endif
          qs   = epsa * es / (txc(i) + epsam1*es)
          tsq  = tx1(i) * tx1(i)
          delq = (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
c
          tx2(i) = tx2(i) + delq
          tx1(i) = tx1(i) - delq * albycp
          tx3(i) = tx3(i) + delq
c
c   third iteration
c
#ifdef ICE
          es   = 10.0*fpvs(tx1(i))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          es   = 10.0*fpvs01
c         es   = 10.0*fpvs0(tx1(i))
#endif
          qs   = epsa * es / (txc(i) + epsam1*es)
          tsq  = tx1(i) * tx1(i)
          delq = (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
c
          tx2(i) = tx2(i) + delq
c
c   "evaporation efficiency" using  kesslers parameterization         *
c
c
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
         ii = i1(i)

         tx2(i) = max(tx2(i) * pii(ii) * dsgbg, 0.0)
         exparg  = zm1p04 * dt * (tx6(i)*tx4(i)*tem) ** zp578
c
         tx3(i) = tx2(i) * (1.0 - exp(exparg)) * min(tx5(i)*rknobt,1.0)
c
         tx3(i) = min(tx6(i), tx3(i))

         tx6(i) = tx6(i) - tx3(i)
#ifdef ISOTOPE
         cond(i)=cond(i)+tmp(i,l)
         condi1(i)=condi1(i)+tmpi1(i,l)
         condi2(i)=condi2(i)+tmpi2(i,l)
         tmp(i,l)=0.D0
         tmpi1(i,l)=0.D0
         tmpi2(i,l)=0.D0
         if (qoi(ii,l).gt.0.D0) then
            virat1=qoi1(ii,l)/qoi(ii,l) !! current vap isotop ratio
            virat2=qoi2(ii,l)/qoi(ii,l)
         else
            virat1=0.D0
            virat2=0.D0
         endif
c         tmp(i,l) = tmp(i,l) - tx3(i) !! tmp should be equal to tx6
#endif
c
         tx3(i)    = tx3(i) * txe(i) * onbdsg(l)
c
         qoi(ii,l) = qoi(ii,l) + tx3(i)
         poi(ii,l) = poi(ii,l) - tx3(i) * pki(ii) * albyp
c
         tcu(ii,l) = tcu(ii,l) - tx3(i) * (albycp*daylen/dt)
         qcu(ii,l) = qcu(ii,l) + tx3(i) * (daylen/dt)
c
#ifdef ISOTOPE
         pl = pii(ii) * sgo(l) * cmb2pa
         tl = poi(ii,l) * pik(ii) * prh(l)
         alpha1=freq(1,tl)
         alpha2=freq(2,tl)
         rh=min(max(f2,1.4*qoi(ii,l)/qs),f1)
c         rh=0.95
c for simple Rayleigh, force rh=0 and nn=0
         if (tl.gt.273.15) then
         gamma1=alpha1*rh/
     $        (alpha1*dif18o**nn*(rh-1)+1.D0)
         beta1=1.D0/(alpha1*dif18o**nn*(1-rh))-1.D0
         gamma2=alpha2*rh/
     $        (alpha2*difhdo**nn*(rh-1)+1)
         beta2=1.D0/(alpha2*difhdo**nn*(1-rh))-1.D0
         else
         gamma1=1.D0
         beta1=0.D0
         gamma2=1.D0
         beta2=0.D0
         endif
c
         flg=0
c
         if (tx3(i).gt.0.D0) then !! evaporation
c
          evap=tx3(i)/txe(i)/onbdsg(l)
c
          tmpi=cond(i)
          cond(i)=max(0.D0,cond(i)-evap)
          evap=max(0.D0,evap-(tmpi-cond(i)))
c
          if (tmpi.gt.0.D0) then
            ff=cond(i)/tmpi
c          
*** Stewart 1975 ***
*   rni/rn=(rni/ro-gamma*qi/q)*ff**beta+gamma*qi/q
*   gamma=alpha_eq*RH/(alpha_eq*(D_iso/D)**n*(RH-1)+1)
*   beta=1/(alpha_eq*(D_iso/D)**n*(1-RH))-1
*   RH: relative humidity
c
            rirat1=condi1(i)/tmpi
            rirat2=condi2(i)/tmpi
            tmpi1o=condi1(i)
            tmpi2o=condi2(i)
c
c lee's equilibrium model
c
            if (eqlee.eq.1) then
              dz=rd*tl/sgo(l)/grav/onbdsg(l)
              call eqm_deg(2,tl,cond(i)/dt*eqlfr,pl,dz,eqfo)
              call eqm_deg(3,tl,cond(i)/dt*eqlfr,pl,dz,eqfd)
c          if ((mype.eq.0).and.(ii.eq.73)) then
c            print'(a6,3i3,6e15.7)',
c     $            'RASEQM',ii,l,flg,tl,cond(i)/dt,pl,dz,eqfo,eqfd
c          endif
c
            else
c constant equilibrium rate
c
              eqfo=eqf
              eqfd=eqf          
            endif
c
            if (l.eq.km-2) then
               eeqfo(ii)=min(eeqfo(ii),eqfo)
               eeqfd(ii)=min(eeqfd(ii),eqfd)
            endif
c     
            condi1(i)=
     $           (rirat1*(eqfo*ff**beta1+1.D0-eqfo)
     $           +eqfo*(1.D0-ff**beta1)*gamma1*virat1)
     $           *cond(i)
            condi2(i)=
     $           (rirat2*(eqfd*ff**beta2+1.D0-eqfd)
     $           +eqfd*(1.D0-ff**beta2)*gamma2*virat2)
     $           *cond(i)
            qoi1(ii,l)=
     $           qoi1(ii,l)+(tmpi1o-condi1(i))*txe(i)*onbdsg(l)
            qoi2(ii,l)=
     $           qoi2(ii,l)+(tmpi2o-condi2(i))*txe(i)*onbdsg(l)
c            
c            if (tl.lt.273.15) then
c               alpha1=1.D0
c               alpha2=1.D0
c            endif
c
c            alp1=(alpha1-1.D0)*eqf+1.D0
c            alp2=(alpha2-1.D0)*eqf+1.D0
c            if (ff.le.f1) then
c               alp1=(alp1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c               alp2=(alp2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c            endif
c
c            qoi1(ii,l)=qoi1(ii,l)+
c     $           condi1(i)*(1.D0-ff**(1./alp1))*txe(i)*onbdsg(l)
c            condi1(i)=condi1(i)*ff**(1./alp1)
c            qoi2(ii,l)=qoi2(ii,l)+
c     $           condi2(i)*(1.D0-ff**(1./alp2))*txe(i)*onbdsg(l)
c            condi2(i)=condi2(i)*ff**(1./alp2)
c
          else
            if ((tmpi.lt.0.D0).or.(condi1(i).ne.0.D0)) then
              print*, i,l,ll,"strange2",tmpi,condi1(i)
            endif
          endif
          if (evap.le.1.D-15) then
             flg=4
             goto 200
          endif
c
          do ll=l+1,k
            cond(i)=cond(i)+tmp(i,ll)
            condi1(i)=condi1(i)+tmpi1(i,ll)
            condi2(i)=condi2(i)+tmpi2(i,ll)
            tmpi=cond(i)
            cond(i)=max(0.D0,cond(i)-evap)
            evap=max(0.D0,evap-(tmpi-cond(i)))
            tmp(i,ll)=0.D0
            tmpi1(i,ll)=0.D0
            tmpi2(i,ll)=0.D0
            if (tmpi.gt.0.D0) then
              ff=cond(i)/tmpi
c
              rirat1=condi1(i)/tmpi
              rirat2=condi2(i)/tmpi
              tmpi1o=condi1(i)
              tmpi2o=condi2(i)
c
              condi1(i)=
     $           (rirat1*(eqfo*ff**beta1+1.D0-eqfo)
     $           +eqfo*(1.D0-ff**beta1)*gamma1*virat1)
     $           *cond(i)
              condi2(i)=
     $           (rirat2*(eqfd*ff**beta2+1.D0-eqfd)
     $           +eqfd*(1.D0-ff**beta2)*gamma2*virat2)
     $           *cond(i)
              qoi1(ii,l)=
     $           qoi1(ii,l)+(tmpi1o-condi1(i))*txe(i)*onbdsg(l)
              qoi2(ii,l)=
     $           qoi2(ii,l)+(tmpi2o-condi2(i))*txe(i)*onbdsg(l)
c            
c              alp1=(alpha1-1.D0)*eqf+1.D0
c              alp2=(alpha2-1.D0)*eqf+1.D0
c              if (ff.le.f1) then
c                 alp1=(alp1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c                 alp2=(alp2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c              endif
c
c              qoi1(ii,l)=qoi1(ii,l)+
c     $             condi1(i)*(1.D0-ff**(1./alp1))*txe(i)*onbdsg(l)
c              condi1(i)=condi1(i)*ff**(1./alp1)
c              qoi2(ii,l)=qoi2(ii,l)+
c     $             condi2(i)*(1.D0-ff**(1./alp2))*txe(i)*onbdsg(l)
c              condi2(i)=condi2(i)*ff**(1./alp2)
c
            else
              if ((tmpi.lt.0.D0).or.(tmpi1(i,ll).ne.0.D0)) then
                print*, i,l,ll,"strange3",tmpi,tmpi1(i,ll)
              endif
            endif
            if (evap.le.1.D-15) then
              flg=4+ll-l
              goto 200
            endif
          enddo
c          print'(a6,3i3,2e15.7)',"EVCN2",mype,ii,l,evap,cond(i)
 200      continue
c
         elseif (tx3(i).lt.0.) then   !! condensation ( q > q' )
c
*** Super saturation Kinetic fractionation ***
* alpha_eff=alpha_eq*alpha_kin
* alpha_kin=S/(alpha_eq*D_iso/D*(S-1)+1)
* S=1-0.003*T [C] (Jouzel and Merlivat 1984)
          if (tl.le.253.15) then
            alpha1=alpha1*
     $            (1-0.003*(tl-273.15))/
     $            (alpha1*dif18o*(-0.003*(tl-273.15))+1)
            alpha2=alpha2*
     $            (1-0.003*(tl-273.15))/
     $            (alpha2*difhdo*(-0.003*(tl-273.15))+1)
*** linear interpolation between -20 to 0 C ***
          elseif (tl.le.273.15) then
            alp0=freq(1,273.15)
            alp20=freq(1,253.15)*
     $           (1-0.003*(-20.))/
     $           (freq(1,253.15)*dif18o*(-0.003*(-20.))+1)
            alpha1=
     $           (alp0*(tl-253.15)
     $           +alp20*(273.15-tl))/20.
            alp0=freq(2,273.15)
            alp20=freq(2,253.15)*
     $           (1-0.003*(-20.))/
     $           (freq(2,253.15)*difhdo*(-0.003*(-20.))+1)
            alpha2=
     $           (alp0*(tl-253.15)
     $           +alp20*(273.15-tl))/20.
          endif
c
         if (qoi(ii,l).ge.0.D0) then !! q > q' >= 0.
         flg=1
c
         if (qoi(ii,l)-tx3(i).gt.0.D0) then
           ff=qoi(ii,l)/(qoi(ii,l)-tx3(i))
         else
           print'(a6,4i3,4e15.7)','wrg4',mype,i,ii,l,
     $           cond(i),cond(i)+tx3(i)/txe(i)/onbdsg(l),
     $           qoi(ii,l),qoi1(ii,l)
           ff=1.D0
         endif
c
c         alpha1=(alpha1-1.D0)*eqf+1.D0
c         alpha2=(alpha2-1.D0)*eqf+1.D0
c         if (ff.le.f1) then
c           alpha1=(alpha1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c           alpha2=(alpha2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c         endif
c
         condi1(i)=condi1(i)+
     $        qoi1(ii,l)*(1.D0-ff**alpha1)
     $        /txe(i)/onbdsg(l)
         cndi1(i)=cndi1(i)+
     $        qoi1(ii,l)*(1.D0-ff**alpha1)
     $        /txe(i)/onbdsg(l)
         qoi1(ii,l)=qoi1(ii,l)*ff**alpha1
         condi2(i)=condi2(i)+
     $        qoi2(ii,l)*(1.D0-ff**alpha2)
     $        /txe(i)/onbdsg(l)
         cndi2(i)=cndi2(i)+
     $        qoi2(ii,l)*(1.D0-ff**alpha2)
     $        /txe(i)/onbdsg(l)
         qoi2(ii,l)=qoi2(ii,l)*ff**alpha2
c     
         elseif (qoi(ii,l).lt.tx3(i)) then !! -q > -q' > 0.
         flg=2
c
         if (qoi(ii,l)-tx3(i).gt.0.D0) then
           ff=qoi(ii,l)/(qoi(ii,l)-tx3(i))
         else
           print'(a6,4i3,4e15.7)','wrg5',mype,i,ii,l,
     $           cond(i),cond(i)+tx3(i)/txe(i)/onbdsg(l),
     $           qoi(ii,l),qoi1(ii,l)
           ff=1.D0
         endif
c
         alpha1=1.D0
         alpha2=1.D0
c         alpha1=(alpha1-1.D0)*eqf+1.D0
c         alpha2=(alpha2-1.D0)*eqf+1.D0
c         if (ff.le.f1) then
c           alpha1=(alpha1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c           alpha2=(alpha2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c         endif
c
         condi1(i)=condi1(i)+
     $        qoi1(ii,l)*(1.D0-ff**alpha1)
     $        /txe(i)/onbdsg(l)
         cndi1(i)=cndi1(i)+
     $        qoi1(ii,l)*(1.D0-ff**alpha1)
     $        /txe(i)/onbdsg(l)
         qoi1(ii,l)=qoi1(ii,l)*ff**alpha1
         condi2(i)=condi2(i)+
     $        qoi2(ii,l)*(1.D0-ff**alpha2)
     $        /txe(i)/onbdsg(l)
         cndi2(i)=cndi2(i)+
     $        qoi2(ii,l)*(1.D0-ff**alpha2)
     $        /txe(i)/onbdsg(l)
         qoi2(ii,l)=qoi2(ii,l)*ff**alpha2
c     
         else !! q >= 0.  q' =< 0
         flg=3
c
         if (cond(i)+qoi(ii,l)/txe(i)/onbdsg(l).gt.0.D0) then
           ff=cond(i)/(cond(i)+qoi(ii,l)/txe(i)/onbdsg(l))
         else
           print'(a6,4i3,4e15.7)','wrg6',mype,i,ii,l,
     $           cond(i),cond(i)+tx3(i)/txe(i)/onbdsg(l),
     $           qoi(ii,l),qoi1(ii,l)
           ff=1.D0
         endif
c
         alpha1=1.D0
         alpha2=1.D0
c         alpha1=(alpha1-1.D0)*eqf+1.D0
c         alpha2=(alpha2-1.D0)*eqf+1.D0
c         if (ff.le.f1) then
c           alpha1=(alpha1-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c           alpha2=(alpha2-1.D0)*(max(ff,f2)-f2)/(f1-f2)+1.D0
c         endif
c
         tmpi=qoi1(ii,l)
         tmpi_2=condi1(i)
         qoi1(ii,l)=
     $        (condi1(i)+tmpi/txe(i)/onbdsg(l))*
     $        (1.D0-ff**(1./alpha1))*txe(i)*onbdsg(l)
         condi1(i)=(condi1(i)+tmpi/txe(i)/onbdsg(l))*ff**(1./alpha1)
         cndi1(i)=cndi1(i)+condi1(i)-tmpi_2
         tmpi=qoi2(ii,l)
         tmpi_2=condi2(i)
         qoi2(ii,l)=
     $        (condi2(i)+tmpi/txe(i)/onbdsg(l))*
     $        (1.D0-ff**(1./alpha2))*txe(i)*onbdsg(l)
         condi2(i)=(condi2(i)+tmpi/txe(i)/onbdsg(l))*ff**(1./alpha2)
         cndi2(i)=cndi2(i)+condi2(i)-tmpi_2
c
         endif
         
         endif
c
*** end for isotope downdraft ***
#endif

      enddo
c
 5000  continue
c
       do i=1,lenb
         ii = i1(i)
         pcu(ii) = pcu(ii) + tx6(i) * (daylen/dt)
         cwmdet(ii) = cwmdet(ii) + cwmw(i)            ! [kg/kg]
#ifdef ISOTOPE
         pcui1(ii) = pcui1(ii) + condi1(i) * (daylen/dt) !! [mm/day]
         pcui2(ii) = pcui2(ii) + condi2(i) * (daylen/dt)
c
         ccnd(ii) = ccnd(ii) + cnd(i) * (daylen/dt)
         ccndi1(ii) = ccndi1(ii) + cndi1(i) * (daylen/dt) !! [mm/day]
         ccndi2(ii) = ccndi2(ii) + cndi2(i) * (daylen/dt)
#endif
       enddo
c
      return
      end
