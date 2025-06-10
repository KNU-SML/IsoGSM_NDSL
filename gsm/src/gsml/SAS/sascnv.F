#define CLIQ cliq
#define CP cp
#define CVAP cvap
#define G g
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define PI pi
#define QMIN qmin
#define RD rd
#define RV rv
#define T0C t0c

      subroutine sascnv(ims2,imx2,kmx,jcap,delt,del,sl,slk,ps,
     &              q1,t1,cldwrk,rn,kbot,ktop,kuo,spd,lat,slimsk,dot,
     &         hpbl,gamt,gamq,qci,qrs,delx)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
c-----------------------------------------------------------------------
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    sascnv      computes convective heating and moisng
c   prgmmr: hua-lu pan       org: w/nmc23    date: 92-03-01
c
c abstract: computes convective heating and moistening using a one
c   cloud type arakawa-schubert convection scheme originally developed
c   by georg grell. the scheme includes updraft and downdraft effects.
c   the closure is the cloud work function. both updraft and downdraft
c   are assumed to be saturated and the heating and moistening are
c   accomplished by the compensating environment. the name comes from
c   "simplified arakawa-schubert convection parameterization".
c
c program history log:
c   92-03-01  hua-lu pan
c
c usage:    call sascnv(im,ix,km,jcap,delt,del,sl,slk,ps,qn,tn,
c    &                q1,t1,rn,kbot,ktop,kuo,spd,lat,slimsk)
c
c   input argument list:
c     im       - integer number of points
c     ix       - leading dimension of qn,tn,q1,t1,spd
c     km       - integer number of levels
c     jcap     - integer spectral truncation
c     dt       - real time step in seconds
c     del      - real (km) sigma layer thickness
c     sl       - real (km) sigma values
c     slk      - real (km) sigma values to the kappa
c     ps       - real (im) surface pressure in kilopascals (cb)
c     qn       - real (ix,km) previous specific humidity in kg/kg
c     tn       - real (ix,km) previous temperature in kelvin
c     q1       - real (ix,km) current specific humidity in kg/kg
c     t1       - real (ix,km) current temperature in kelvin
c     spd      - real (ix,km) current wind speed
c     lat      - integer  current latitude index
c     slimsk   - real (im) land(1),sea(0), ice(2) flag
c
c   output argument list:
c     q1       - real (ix,km) adjusted specific humidity in kg/kg
c     t1       - real (ix,km) adjusted temperature in kelvin
c     rn       - real (im) convective rain in meters
c     kbot     - integer (im) cloud bottom level
c     ktop     - integer (im) cloud top level
c     kuo      - integer (im) bit flag indicating deep convection
c
c subprograms called:
c   fpvs     - function to compute saturation vapor pressure
c
c remarks: function fpvs is inlined by fpp.
c          nonstandard automatic arrays are used.
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
c$$$
#include <paramodel.h>
      dimension del(imx2,kmx),sl(imx2,kmx),slk(imx2,kmx),ps(imx2),
c     &          q1(imx2,kmx*ncld),t1(imx2,kmx),rn(imx2),
     &          q1(imx2,kmx),t1(imx2,kmx),rn(imx2),
     &          kbot(imx2),ktop(imx2),kuo(imx2)
      dimension spd(imx2,kmx),slimsk(imx2),dot(imx2,kmx),cldwrk(imx2)
      dimension hpbl(imx2),gamt(imx2),gamq(imx2)
      dimension qci(imx2,kmx),qrs(imx2,kmx)
c  physical parameters
      parameter(cp=cp_,g=g_,rd=rd_,rv=rv_,hvap=hvap_,t0c=t0c_)
      parameter(cvap=cvap_,cliq=cliq_)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,
     &          el2orc=hvap*hvap/(rv*cp),eps=rd/rv,epsm1=rd/rv-1.,
     &          fv=1./eps-1.)
      parameter(terr=0.,c0=.002,qmin=qmin_)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
c
c  local variables and arrays
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      dimension p(ILOTS,levs_),pdot(ILOTS),
     *          acrtfct(ILOTS),
     &          to(ILOTS,levs_),
     *          qo(ILOTS,levs_)
      dimension qeso(ILOTS,levs_)
      dimension tvo(ILOTS,levs_),
     *          dbyo(ILOTS,levs_),
     &          zo(ILOTS,levs_),
     &          heo(ILOTS,levs_),
     *          heso(ILOTS,levs_),
     &          qrcd(ILOTS,levs_),
     &          dellah(ILOTS,levs_),
     *          dellaq(ILOTS,levs_),
     &          hcko(ILOTS,levs_),
     &          qcko(ILOTS,levs_),
     *          eta(ILOTS,levs_),
     &          etad(ILOTS,levs_),
     &          qrcdo(ILOTS,levs_),
     &          pwo(ILOTS,levs_),
     *          pwdo(ILOTS,levs_),
     *          dtconv(ILOTS),
     &          deltv(ILOTS),acrt(ILOTS)
      dimension psfc(ILOTS),hmax(ILOTS),
     *          kb(ILOTS), delq(ILOTS),
     &          hkbo(ILOTS),qkbo(ILOTS),
     *          kbcon(ILOTS),pbcdif(ILOTS),
     &          vmax(ILOTS),kds(ILOTS),
     &          hmin(ILOTS),lmin(ILOTS),
     *          jmin(ILOTS),pwavo(ILOTS),
     &          aa1(ILOTS),vshear(ILOTS),
     *          shrmax(ILOTS),
     &          kshmax(ILOTS),edt(ILOTS),
     &          edto(ILOTS),pwevo(ILOTS),
     &          qcond(ILOTS),
     &          hcdo(ILOTS),qcdo(ILOTS),
     *          ddp(ILOTS),pp2(ILOTS),
     &          adet(ILOTS),aatmp(ILOTS),
     &          xhkb(ILOTS),xqkb(ILOTS),
     &          xpwav(ILOTS),xpwev(ILOTS),
     *          xhcd(ILOTS),
     &          xaa0(ILOTS),f(ILOTS),
     *          xk(ILOTS),
     &          xmb(ILOTS),ktcon(ILOTS),
     &          edtx(ILOTS),xqcd(ILOTS),
     &          hsbar(ILOTS),xmbmax(ILOTS),
     &          xlamb(ILOTS),xlamd(ILOTS),
     &          kbdtr(ILOTS),
     &          excess(ILOTS),kpbl(ILOTS),
     &          plcl(ILOTS),klcl(ILOTS),
     *          ktdown(ILOTS)
      dimension delhbar(ILOTS),delqbar(ILOTS),
     *          deltbar(ILOTS)
      dimension pcrit(15), acritt(15), acrit(15)
c
      real mbdt
      save pcrit, acritt
      logical totflg,
     *        cnvflg(ILOTS),dwnflg(ILOTS),
     1        dwnflg2(ILOTS), flg(ILOTS),
     &        lfcflg,gcmflg,lclflg,evaflg
      data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,
     &           350.,300.,250.,200.,150./
      data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
ccc   parameter(xk1=.2,xlhor=10000.,xhver=5000.,
      parameter(xk1=2.e-5,xlhor=3.e4,xhver=5000.,theimax=1.,
ccc  &          xc1=5.e-7,xc2=3.e4,xc3=10.,ecesscr=3.0)
     &          xc1=1.e-7,xc2=1.e4,xc3=3.e3,ecesscr=3.0,edtk1=3.e4)
c  gdas derived acrit
c     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
c    &            .743,.813,.886,.947,1.138,1.377,1.896/
c-----------------------------------------------------------------------
c  initialize arrays
c
      im=ims2
      km=kmx
c
#ifdef RSM
        gcmflg = .false.
        betal = .3
        evaflg = .false.   !!! evaporational efficiency
#else
        gcmflg = .true.
        betal = .05
        evaflg = .true.   !!! evaporational efficiency
#endif
c
      lfcflg = .true.   !!! old trigger
      lclflg = .true.   !!! lcl cloud bottom in the new trigger
c
c     ncloud = ncld
      ncloud = 1
      kci = 0
      krs = 0
      if(ncloud.gt.1) then
         kci = km
         krs = 2*km
      endif
      do i=1,im
        rn(i)=0.
        kbot(i)=km+1
        ktop(i)=0
        kuo(i)=0
        cnvflg(i) = .true.
        dtconv(i) = 3600.
        cldwrk(i) = 0.
        pdot(i) = 0.
        xmbmax(i) = .3
        excess(i) = 0.0
        plcl(i) = 0.0
        kpbl(i) = 1
      enddo
      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo
      dt2 = 2. * delt
      dtmin = max(dt2,1200.)
      dtmax = max(dt2,3600.)
c  model tunable parameters are all here
      mbdt = 10.
      if(gcmflg) then
         edtmaxl = min(edtk1/delx,.9)
         edtmaxs = .9
         if(evaflg) edtmaxs = .3
      else
         edtmaxl = .3
         edtmaxs = .3
      endif
      alphal = .5
      alphas = .5
      betas = .05
      evef = 0.07
      if(evaflg) then
         evfact = 1. - edtmaxs
      endif
      pdpdwn=0.
      pdetrn=200.
      w1l = -2.e-3 * (jcap / 62.)
      w2l = -1.e-2 * (jcap / 62.)
ccc      w3l = -1.e-2 * (jcap / 62.)
ccc      w4l = -1.e-3 * (jcap / 62.)
      w3l = -2.e-3 * (jcap / 62.)
      w4l = -2.e-4 * (jcap / 62.)
      w1s = -2.e-3 * (jcap / 62.)
      w2s = -1.e-2 * (jcap / 62.)
      w3s = -2.e-3 * (jcap / 62.)
      w4s = -2.e-4 * (jcap / 62.)
ccccc if(im.eq.384) then
ccccc   latd = 45
ccccc   lond = 376
ccccc elseif(im.eq.768) then
ccccc   latd = 80
ccccc   lond = 81
ccccc else
ccccc   latd = 0
ccccc   lond = 0
ccccc endif
c
c  define top layer for search of the downdraft originating layer
c  and the maximum thetae for updraft
c
c      kbmax = km
c      kbm = km
c      kmax = km
      kbmax = 1
      kbm = 1
      kmax = 1
      do k = 1, km
        do i = 1, im
        if((sl(i,k).gt..45).and.(k.gt.kbmax-1)) kbmax = k + 1
        if((sl(i,k).gt..7).and.(k.gt.kbm-1)) kbm = k + 1
        if((sl(i,k).gt..05).and.(k.gt.kmax-1)) kmax = k + 1
        enddo
      enddo
      kmax = min(kmax,km)
c
c   convert surface pressure to mb from cb
c
      do i = 1, im
        psfc(i) = ps(i) * 10.
      enddo
      do k = 1, km
        do i = 1, im
          qci(i,k) = 0.0
          qrs(i,k) = 0.0
          pwo(i,k) = 0.
          pwdo(i,k) = 0.
        enddo
      enddo
      do k = 1, kmax
        do i = 1, im
          p(i,k) = psfc(i) * sl(i,k)
          pwo(i,k) = 0.
          pwdo(i,k) = 0.
          to(i,k) = t1(i,k)
          qo(i,k) = q1(i,k)
          dbyo(i,k) = 0.
        enddo
      enddo
c
c  column variables
c  p is pressure of the layer (mb)
c  t is temperature at t-dt (k)..tn
c  q is mixing ratio at t-dt (kg/kg)..qn
c  to is temperature at t+dt (k)... this is after advection and turbulan
c  qo is mixing ratio at t+dt (kg/kg)..q1
c
      do k = 1, kmax
        do i = 1, im
#ifdef ICE
          qeso(i,k) = 10. * fpvs(t1(i,k))
#else
          qeso(i,k) = 10. * fpvs0(t1(i,k))
#endif
          qeso(i,k) = eps * qeso(i,k) / (p(i,k) + epsm1 * qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
          tvo(i,k) = to(i,k) + fv * to(i,k) * max(qo(i,k),qmin)
        enddo
      enddo
c
c  hydrostatic height assume zero terr
c
      do i = 1, im
        dlnsig = log(sl(i,1))
        zo(i,1) = terr - dlnsig * rd / g * tvo(i,1)
      enddo
      do k = 2, kmax
        do i = 1, im
          dlnsig = log(sl(i,k) / sl(i,k-1))
          zo(i,k) = zo(i,k-1) - dlnsig * rd / g
     &              * .5 * (tvo(i,k) + tvo(i,k-1))
        enddo
      enddo
c  compute moist static energy
      do k = 1, kmax
        do i = 1, im
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
c         heo(i,k) = min(heo(i,k),heso(i,k))
        enddo
      enddo
c
c  determine level with largest moist static energy
c  this is the level where updraft starts
c
      do i = 1, im
        hmax(i) = heo(i,1)
        kb(i) = 1
      enddo
      do k = 2, kbm
        do i = 1, im
          if(heo(i,k).gt.hmax(i).and.cnvflg(i)) then
            kb(i) = k
            hmax(i) = heo(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(qo(i,kb(i)).lt.qmin) cnvflg(i) = .false.
      enddo
c
c  search for downdraft originating level above theta-e minimum
c
      do i = 1, im
        hmin(i) = heso(i,1)
        lmin(i) = kbmax
        jmin(i) = kbmax
      enddo
      do k = 2, kbmax
        do i = 1, im
          if(heso(i,k).lt.hmin(i).and.cnvflg(i)) then
             lmin(i) = k + 1
             hmin(i) = heso(i,k)
          endif
        enddo
      enddo
      do k = 1, kmax - 1
        do i = 1, im
        if(cnvflg(i)) then
          dz = .5 * (zo(i,k+1) - zo(i,k))
          dp = .5 * (p(i,k+1) - p(i,k))
#ifdef ICE
          es = 10. * fpvs(to(i,k+1))
#else
          es = 10. * fpvs0(to(i,k+1))
#endif
          pprime = p(i,k+1) + epsm1 * es
          qs = eps * es / pprime
          dqsdp = - qs / pprime
          desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
          dqsdt = qs * p(i,k+1) * desdt / (es * pprime)
          gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
          dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
          dq = dqsdt * dt + dqsdp * dp
          to(i,k) = to(i,k+1) + dt
          qo(i,k) = qo(i,k+1) + dq
          po = .5 * (p(i,k) + p(i,k+1))
#ifdef ICE
          qeso(i,k) = 10. * fpvs(to(i,k))
#else
          qeso(i,k) = 10. * fpvs0(to(i,k))
#endif
          qeso(i,k) = eps * qeso(i,k) / (po + epsm1 * qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
c         qo(i,k) = min(qo(i,k),qeso(i,k))
          heo(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                cp * to(i,k) + hvap * qeso(i,k)
        endif
        enddo
      enddo
c
c  look for convective cloud base as the level of free convection
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hkbo(i) = heo(i,indx)
          qkbo(i) = qo(i,indx)
        endif
      enddo
      do i = 1, im
         flg(i) = cnvflg(i)
         kbcon(i) = kmax
      enddo
c
      if(lfcflg) then
          do k = 1, kbmax
            do i = 1, im
              if(flg(i).and.k.gt.kb(i)) then
                hsbar(i) = heso(i,k)
                if(hkbo(i).gt.hsbar(i)) then
                  flg(i) = .false.
                  kbcon(i) = k
                endif
              endif
            enddo
          enddo
          do i = 1, im
            if(cnvflg(i)) then
              pbcdif(i) = -p(i,kbcon(i)) + p(i,kb(i))
c             pdot(i) = 10.* dot(i,kbcon(i))
              if(pbcdif(i).gt.150.) cnvflg(i) = .false.
              if(kbcon(i).eq.kmax) cnvflg(i) = .false.
            endif
          enddo
      else
c
c  determine mesoscale trigger
c
      do k = 1,kmax
         do i = 1,im
            if(flg(i).and.zo(i,k).ge.hpbl(i)) then
               kpbl(i) = k
               flg(i) = .false.
            endif
         enddo
      enddo
      do i = 1,im
         if(cnvflg(i)) then
            indx = kb(i)
            indp = kpbl(i)
            thei = theimax/pi_*(atan(xk1*(delx-xlhor))+pi_/2.)
            thec = gamt(i) + eps * tvo(i,1) * gamq(i)
            if(zo(i,indx).le.hpbl(i)) then
               theb = thec
            else
               delza = max(zo(i,indx) - hpbl(i),qmin)
               thekb = tvo(i,indx)/slk(i,indx)
               thekh = tvo(i,indp)/slk(i,indp)
               theavg = (thekb+thekh)/2.
               thedif = (thekb-thekh)/delza
               dzfac = min(xc1*delza**2.+xc2/theavg*thedif,50.)
               dzfac = max(dzfac,0.)
               theb = thec*exp(-dzfac)
            endif
            omgkb = dot(i,indx)*10.
            omgkbp1 = dot(i,indx+1)*10.
            omgdif = (omgkb-omgkbp1)/(p(i,indx)-p(i,indx+1))
            omgfac = (min(abs(omgdif),1.e50))**(1./3.)
            if(omgdif.lt.0.) then
               omgfac = -1.*omgfac
            endif
            theom = max((thei+theb)*(1.+delx/xc3*omgfac),0.)
            excess(i) = min(theom,ecesscr)
         endif
      enddo
c
      do i = 1,im
         if(cnvflg(i)) then
            indx = kb(i)
            rh = min(max(qo(i,indx)/qeso(i,indx),qmin),1.)
            thermal = to(i,indx)+excess(i)
            chi = thermal/(1669.0-122.0*rh-thermal)
            plcl(i) = p(i,indx)*(rh**chi)
            if(plcl(i).lt.p(i,kbmax).or.plcl(i).gt.p(i,indx)) then
               cnvflg(i) = .false.
            endif
         endif
      enddo
      do i = 1, im
         flg(i) = cnvflg(i)
         klcl(i) = kmax
      enddo
      do k = 1,kbmax
         do i = 1,im
            if(flg(i).and.k.gt.kb(i).and.p(i,k).le.plcl(i)) then
               klcl(i) = k
               flg(i) = .false.
            endif
         enddo
      enddo
c
      do i = 1,im
         if(cnvflg(i)) then
            if(klcl(i).ge.kbmax) cnvflg(i) = .false.
            if((hkbo(i)+excess(i)*cp).lt.heso(i,klcl(i)))
     1          cnvflg(i) = .false.
         endif
      enddo
c
c     upto here, kbcon is lcl. put the cloud bottom as lfc
c
      do i = 1,im
         if((.not.lclflg).and.cnvflg(i)) then
            flg(i) = .true.
         else
            kbcon(i) = klcl(i)
         endif
      enddo
c
      do k = 1, kbmax
         do i = 1, im
            if((.not.lclflg).and.cnvflg(i)) then
               if(flg(i).and.k.gt.kb(i).and.k.ge.klcl(i)) then
                  hsbar(i) = heso(i,k)
                  if((hkbo(i)+excess(i)*cp).gt.hsbar(i)) then
                     flg(i) = .false.
                     kbcon(i) = k
                  endif
               endif
            endif
         enddo
      enddo
c
      endif
c
      do i = 1,im
         if(cnvflg(i)) then
            if(kbcon(i).ge.kbmax) cnvflg(i) = .false.
            pdot(i) = 10.* dot(i,kbcon(i))
         endif
      enddo
c
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c  found lfc, can define rest of variables
 6001 format(2x,-2p10f12.2)
 6002 format(2x,10f12.2)
 6003 format(2x,3p10f12.2)
c
c  determine entrainment rate between kb and kbcon
c
      do i = 1, im
        alpha = alphas
        if(slimsk(i).eq.1.) alpha = alphal
        if(cnvflg(i)) then
          if(kb(i).eq.1) then
            dz = .5 * (zo(i,kbcon(i)) + zo(i,kbcon(i)-1)) - zo(i,1)
          else
            dz = .5 * (zo(i,kbcon(i)) + zo(i,kbcon(i)-1))
     &         - .5 * (zo(i,kb(i)) + zo(i,kb(i)-1))
          endif
          if(kbcon(i).ne.kb(i)) then
            xlamb(i) = -log(alpha) / dz
          else
            xlamb(i) = 0.
          endif
        endif
      enddo
c  determine updraft mass flux
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            eta(i,k) = 1.
          endif
        enddo
      enddo
      do k = kbmax, 2, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.kbcon(i).and.k.ge.kb(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            eta(i,k) = eta(i,k+1) * exp(-xlamb(i) * dz)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.kb(i).eq.1.and.kbcon(i).gt.1) then
          dz = .5 * (zo(i,2) - zo(i,1))
          eta(i,1) = eta(i,2) * exp(-xlamb(i) * dz)
        endif
      enddo
c
c  work up updraft cloud properties
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = hkbo(i)
          qcko(i,indx) = qkbo(i)
          pwavo(i) = 0.
        endif
      enddo
c
c  cloud property below cloud base is modified by the entrainment proces
c
      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.le.kbcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
          if(cnvflg(i).and.k.gt.kbcon(i)) then
            hcko(i,k) = hcko(i,k-1)
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
        enddo
      enddo
c  determine cloud top
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon(i) = 1
      enddo
c     do k = 2, kmax
c       kk = kmax - k + 1
c       do i = 1, im
c         if(dbyo(i,kk).ge.0..and.flg(i).and.kk.gt.kbcon(i)) then
c           ktcon(i) = kk + 1
c           flg(i) = .false.
c         endif
c       enddo
c     enddo
c
c     check inversion
c
      do k = 2, kmax
        do i = 1, im
ccc       if(dbyo(i,k).lt.0..and.flg(i).and.k.gt.
          if((dbyo(i,k)+excess(i)*cp).lt.0..and.flg(i).and.k.gt.
     1        kbcon(i)) then
            ktcon(i) = k
            flg(i) = .false.
          endif
        enddo
      enddo
      do i = 1, im
        flg(i) = cnvflg(i)
        ktdown(i) = 1
      enddo
      do k = kmax,1,-1
        do i = 1, im
          if(dbyo(i,k).ge.0..and.flg(i).and.k.gt.
     1        kbcon(i).and.k.le.ktcon(i)) then
            ktdown(i) = k
            flg(i) = .false.
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
            ktcon(i) = ktdown(i) + 1
        endif
      enddo
c
c    check cloud depth
c
      do i = 1, im
        if(cnvflg(i).and.(p(i,kbcon(i)) - p(i,ktcon(i))).lt.150.)
     &  cnvflg(i) = .false.
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  detraining cloud
c
      do i = 1, im
        if(cnvflg(i)) then
          dz = .5 * (zo(i,ktcon(i))+zo(i,ktcon(i)-1))
     &       - .5 * (zo(i,kbcon(i))+zo(i,kbcon(i)-1))
          xlamb(i) = -1. / dz
        endif
        dwnflg(i) = cnvflg(i)
        if(cnvflg(i).and.(p(i,kbcon(i))-p(i,ktcon(i))).gt.pdetrn)
     &     dwnflg(i)=.false.
        dwnflg2(i) = cnvflg(i)
        if(cnvflg(i).and.(p(i,kbcon(i))-p(i,ktcon(i))).lt.pdpdwn)
     &     dwnflg2(i)=.false.
      enddo
      do k = 2, kmax - 1
        do i = 1, im
          if(dwnflg(i).and.k.gt.kbcon(i).and.k.le.ktcon(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            eta(i,k) = eta(i,k-1) * exp( xlamb(i) * dz)
          endif
        enddo
      enddo
c
c  cloud property above cloud top is modified by the detrainment process
c
      do k = 2, kmax - 1
        do i = 1, im
          if(dwnflg(i).and.k.gt.kbcon(i).and.k.le.ktcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
            dbyo(i,k) = hcko(i,k) - heso(i,k)
          endif
        enddo
      enddo
c
c  make sure that jmin is within the cloud
c
      do i = 1, im
        if(cnvflg(i)) then
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon(i)+1)
        endif
      enddo
      do i = 1, im
        if(cnvflg(i).and.dwnflg2(i).and.jmin(i).le.kbcon(i))
     &     then
          cnvflg(i) = .false.
          dwnflg(i) = .false.
          dwnflg2(i) = .false.
        endif
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  compute cloud moisture property and precipitation
c
      do i = 1, im
          aa1(i) = 0.
      enddo
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.lt.ktcon(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            dz1 = (zo(i,k) - zo(i,k-1))
            gamma = el2orc * qeso(i,k) / (to(i,k)**2)
            qrch = qeso(i,k)
     &           + gamma * dbyo(i,k) / (hvap * (1. + gamma))
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            qcko(i,k) = factor * qcko(i,k-1) + onemf *
     &                  .5 * (qo(i,k) + qo(i,k+1))
            dq = eta(i,k) * qcko(i,k) - eta(i,k) * qrch
c
c  below lfc check if there is excess moisture to release latent heat
c
            if(dq.gt.0.) then
              etah = .5 * (eta(i,k) + eta(i,k-1))
              qlk = dq / (eta(i,k) + etah * c0 * dz)
              aa1(i) = aa1(i) - dz1 * g * qlk
              qc = qlk + qrch
              pwo(i,k) = etah * c0 * dz * qlk
              qcko(i,k) = qc
              pwavo(i) = pwavo(i) + pwo(i,k)
              qci(i,k) = qlk
            endif
          endif
        enddo
      enddo
c
c  calculate cloud work function at t+dt
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.kbcon(i).and.k.lt.ktcon(i)) then
            dz1 = zo(i,k) - zo(i,k-1)
            gamma = el2orc * qeso(i,k-1) / (to(i,k-1)**2)
            rfact =  1. + fv * cp * gamma
     &               * to(i,k-1) / hvap
            aa1(i) = aa1(i) +
     &               dz1 * (g / (cp * to(i,k-1)))
     &               * dbyo(i,k-1) / (1. + gamma)
     &               * rfact
            aa1(i)=aa1(i)+
     &               dz1 * g * fv *
     &               max(0.,(qeso(i,k-1) - qo(i,k-1)))
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) dwnflg(i) = .false.
        if(cnvflg(i).and.aa1(i).le.0.) dwnflg2(i) = .false.
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
c
c------- downdraft calculations
c
c  determine level with largest wind speed
      do i = 1, im
          vmax(i) = spd(i,kb(i))
          kds(i) = kb(i)
      enddo
      do k = 2, kmax - 1
        do i = 1, im
          if(k.ge.kb(i).and.spd(i,k).gt.vmax(i).and.cnvflg(i)
     &       .and.k.le.ktcon(i)) then
            vmax(i) = spd(i,k)
            kds(i) = k
          endif
        enddo
      enddo
c
c--- determine downdraft strength in terms of windshear
c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
          shrmax(i) = 0.
          kshmax(i) = 1
        endif
      enddo
      do k = 1, kmax
        do i = 1, im
cgcm          if(k.ge.kb(i).and.k.le.kds(i).and.cnvflg(i)) then
          if(k.ge.kb(i).and.k.le.ktcon(i).and.cnvflg(i)) then
            shear = abs((spd(i,k+1)-spd(i,k))/(zo(i,k+1)-zo(i,k)))
            vshear(i) = vshear(i) + shear
            shrmax(i) = max(shrmax(i),shear)
            if(shrmax(i).eq.shear) kshmax(i) = k
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
cgcm          knumb = kds(i) - kb(i) + 1
          knumb = ktcon(i) - kb(i) + 1
          knumb = max(knumb,1)
          vshear(i) = 1.e3 * vshear(i) / float(knumb)
          if(vshear(i).ge.1.35) then
          e1 = 1.591-.639*vshear(i)
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          else
          e1 = 0.9
          endif
          edt(i)=1.-e1
          edt(i) = min(edt(i),1.0)
          edt(i) = max(edt(i),.1)
          edto(i)=edt(i)
          edtx(i)=edt(i)
        endif
      enddo
c  determine detrainment rate between 1 and kbdtr
      do i = 1, im
        kbdtr(i) = kbcon(i)
        beta = betas
        if(slimsk(i).eq.1.) beta = betal
        if(cnvflg(i)) then
          kbdtr(i) = kbcon(i)
          kbdtr(i) = max(kbdtr(i),1)
          xlamd(i) = 0.
          if(kbdtr(i).gt.1) then
            dz = .5 * zo(i,kbdtr(i)) + .5 * zo(i,kbdtr(i)-1)
     &         - zo(i,1)
            xlamd(i) = log(beta) / dz
          endif
        endif
      enddo
c  determine downdraft mass flux
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            etad(i,k) = 1.
          endif
          qrcdo(i,k) = 0.
        enddo
      enddo
      do k = kbmax, 2, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.kbdtr(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            etad(i,k) = etad(i,k+1) * exp(xlamd(i) * dz)
          endif
        enddo
      enddo
      k = 1
      do i = 1, im
        if(cnvflg(i).and.kbdtr(i).gt.1) then
          dz = .5 * (zo(i,2) - zo(i,1))
          etad(i,k) = etad(i,k+1) * exp(xlamd(i) * dz)
        endif
      enddo
c
c--- downdraft moisture properties
c
      do i = 1, im
        pwevo(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i) = heo(i,jmn)
          qcdo(i) = qo(i,jmn)
          qrcdo(i,jmn) = qeso(i,jmn)
        endif
      enddo
      do k = kmax-1, 1, -1
        do i = 1, im
          if(cnvflg(i).and.k.lt.jmin(i)) then
            dq = qeso(i,k)
            dt = to(i,k)
            gamma = el2orc * dq / dt**2
            dh = hcdo(i) - heso(i,k)
            qrcdo(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
            detad = etad(i,k+1) - etad(i,k)
            pwdo(i,k) = etad(i,k+1) * qcdo(i) -
     &                 etad(i,k) * qrcdo(i,k)
              pwdo(i,k) = pwdo(i,k) - detad *
     &                   .5 * (qrcdo(i,k) + qrcdo(i,k+1))
            qcdo(i) = qrcdo(i,k)
            pwevo(i) = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo
c
c--- final downdraft strength dependent on precip
c--- efficiency (edt), normalized condensate (pwav), and
c--- evaporate (pwev)
c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(dwnflg2(i)) then
          if(pwevo(i).lt.0.) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            if(gcmflg.and.edto(i).ge..90) then
               cnvflg(i) = .false.
               dwnflg2(i) = .false.
            endif
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = 0.
          endif
        else
          edto(i) = 0.
        endif
      enddo
c
c
c--- downdraft cloudwork functions
c
c
      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            gamma = el2orc * qeso(i,k) / to(i,k)**2
            dhh=hcdo(i)
            dt=to(i,k)
            dg=gamma
            dh=heso(i,k)
            dz=-1.*(zo(i,k+1)-zo(i,k))
            aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))
     &             *(1.+fv*cp*dg*dt/hvap)
            aa1(i)=aa1(i)+edto(i)*
     &      dz*g*fv*max(0.,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
      do i = 1, im
        if(aa1(i).le.0.) cnvflg(i) = .false.
        if(aa1(i).le.0.) dwnflg(i) = .false.
        if(aa1(i).le.0.) dwnflg2(i) = .false.
      enddo
c
c
c--- what would the change be, that a cloud with unit mass
c--- will do to the environment?
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          dp = 100. * psfc(i) * del(i,1)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i)
     &                - heo(i,1)) * g / dp
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i)
     &                - qo(i,1)) * g / dp
        endif
      enddo
c
c--- changed due to subsidence and entrainment
c
      do k = 2, kmax-1
        do i = 1, im
          if(cnvflg(i).and.k.lt.ktcon(i)) then
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            dv1= heo(i,k)
            dv2 = .5 * (heo(i,k) + heo(i,k+1))
            dv3= heo(i,k-1)
            dv1q= qo(i,k)
            dv2q = .5 * (qo(i,k) + qo(i,k+1))
            dv3q= qo(i,k-1)
            dp = 100. * psfc(i) * del(i,k)
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            deta = eta(i,k) - eta(i,k-1)
            detad = etad(i,k) - etad(i,k-1)
            dellah(i,k) = dellah(i,k) +
     &          ((aup * eta(i,k) - adw * edto(i) * etad(i,k)) * dv1
     &      - (aup * eta(i,k-1) - adw * edto(i) * etad(i,k-1))* dv3
     &                   - aup * deta * dv2
     &                  + adw * edto(i) * detad * hcdo(i)) * g / dp
            dellaq(i,k) = dellaq(i,k) +
     &          ((aup * eta(i,k) - adw * edto(i) * etad(i,k)) * dv1q
     &      - (aup * eta(i,k-1) - adw * edto(i) * etad(i,k-1))* dv3q
     &                  - aup * deta * dv2q
     &     +adw*edto(i)*detad*.5*(qrcdo(i,k)+qrcdo(i,k-1))) * g / dp
          endif
        enddo
      enddo
c
c------- cloud top
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 100. * psfc(i) * del(i,indx)
          dv1 = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *
     &                     (hcko(i,indx-1) - dv1) * g / dp
          dvq1 = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *
     &                     (qcko(i,indx-1) - dvq1) * g / dp
        endif
      enddo
c
c------- final changed variable per unit mass flux
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.gt.ktcon(i)) then
            qo(i,k) = q1(i,k)
            to(i,k) = t1(i,k)
          endif
          if(cnvflg(i).and.k.le.ktcon(i)) then
            qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
            dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
            to(i,k) = dellat * mbdt + t1(i,k)
            qo(i,k) = max(qo(i,k),qmin)
          endif
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c--- the above changed environment is now used to calulate the
c--- effect the arbitrary cloud (with unit mass flux)
c--- would have on the stability,
c--- which then is used to calculate the real mass flux,
c--- necessary to keep this change in balance with the large-scale
c--- destabilization.
c
c--- environmental conditions again, first heights
c
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i)) then
#ifdef ICE
            qeso(i,k) = 10. * fpvs(to(i,k))
#else
            qeso(i,k) = 10. * fpvs0(to(i,k))
#endif
            qeso(i,k) = eps * qeso(i,k) / (p(i,k) + epsm1 * qeso(i,k))
            qeso(i,k) = max(qeso(i,k),qmin)
            tvo(i,k) = to(i,k) + fv * to(i,k) * max(qo(i,k),qmin)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
          xpwav(i) = 0.
        endif
      enddo
c
c  hydrostatic height assume zero terr
c
      do i = 1, im
        dlnsig = log(sl(i,1))
        if(cnvflg(i)) then
          zo(i,1) = terr - dlnsig * rd / g * tvo(i,1)
        endif
      enddo
      do k = 2, kmax
        do i = 1, im
          dlnsig = log(sl(i,k) / sl(i,k-1))
          if(cnvflg(i)) then
            zo(i,k) = zo(i,k-1) - dlnsig * rd / g
     &             * .5 * (tvo(i,k) + tvo(i,k-1))
          endif
        enddo
      enddo
c
c--- moist static energy
c
      do k = 1, kmax - 1
        do i = 1, im
          if(cnvflg(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (p(i,k+1) - p(i,k))
#ifdef ICE
            es = 10. * fpvs(to(i,k+1))
#else
            es = 10. * fpvs0(to(i,k+1))
#endif
            pprime = p(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * p(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po = .5 * (p(i,k) + p(i,k+1))
#ifdef ICE
            qeso(i,k) = 10. * fpvs(to(i,k))
#else
            qeso(i,k) = 10. * fpvs0(to(i,k))
#endif
            qeso(i,k) = eps * qeso(i,k) / (po + epsm1 * qeso(i,k))
            qeso(i,k) = max(qeso(i,k),qmin)
c           qo(i,k) = min(qo(i,k),qeso(i,k))
            heo(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
      k = kmax
      do i = 1, im
        if(cnvflg(i)) then
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
c         heo(i,k) = min(heo(i,k),heso(i,k))
        endif
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          xhkb(i) = heo(i,indx)
          xqkb(i) = qo(i,indx)
          hcko(i,indx) = xhkb(i)
          qcko(i,indx) = xqkb(i)
        endif
      enddo
c
c**************************** static control
c
c------- moisture and cloud work functions
c
      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.le.ktcon(i)) then
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            hcko(i,k) = factor * hcko(i,k-1) + onemf *
     &                  .5 * (heo(i,k) + heo(i,k+1))
          endif
        enddo
      enddo
      do k = 2, kmax - 1
        do i = 1, im
          if(cnvflg(i).and.k.gt.kb(i).and.k.lt.ktcon(i)) then
            dz = .5 * (zo(i,k+1) - zo(i,k-1))
            gamma = el2orc * qeso(i,k) / (to(i,k)**2)
            xdby = hcko(i,k) - heso(i,k)
            xdby = max(xdby,0.)
            xqrch = qeso(i,k)
     &           + gamma * xdby / (hvap * (1. + gamma))
            factor = eta(i,k-1) / eta(i,k)
            onemf = 1. - factor
            qcko(i,k) = factor * qcko(i,k-1) + onemf *
     &                  .5 * (qo(i,k) + qo(i,k+1))
            dq = eta(i,k) * qcko(i,k) - eta(i,k) * xqrch
            if(dq.gt.0.) then
              etah = .5 * (eta(i,k) + eta(i,k-1))
              qlk = dq / (eta(i,k) + etah * c0 * dz)
              xaa0(i) = xaa0(i) - (zo(i,k) - zo(i,k-1)) * g * qlk
              xqc = qlk + xqrch
              xpw = etah * c0 * dz * qlk
              qcko(i,k) = xqc
              xpwav(i) = xpwav(i) + xpw
            endif
          endif
          if(cnvflg(i).and.k.gt.kbcon(i).and.k.lt.ktcon(i)) then
            dz1 = zo(i,k) - zo(i,k-1)
            gamma = el2orc * qeso(i,k-1) / (to(i,k-1)**2)
            rfact =  1. + fv * cp * gamma
     &               * to(i,k-1) / hvap
            xdby = hcko(i,k-1) - heso(i,k-1)
            xaa0(i) = xaa0(i)
     &              + dz1 * (g / (cp * to(i,k-1)))
     &              * xdby / (1. + gamma)
     &              * rfact
            xaa0(i)=xaa0(i)+
     &               dz1 * g * fv *
     &               max(0.,(qeso(i,k-1) - qo(i,k-1)))
          endif
        enddo
      enddo
c
c------- downdraft calculations
c
c
c--- downdraft moisture properties
c
      do i = 1, im
        xpwev(i) = 0.
      enddo
      do i = 1, im
        if(dwnflg2(i)) then
          jmn = jmin(i)
          xhcd(i) = heo(i,jmn)
          xqcd(i) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
        endif
      enddo
      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            dq = qeso(i,k)
            dt = to(i,k)
            gamma = el2orc * dq / dt**2
            dh = xhcd(i) - heso(i,k)
            qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
            detad = etad(i,k+1) - etad(i,k)
            xpwd = etad(i,k+1) * qrcd(i,k+1) -
     &                 etad(i,k) * qrcd(i,k)
            xpwd = xpwd - detad *
     &             .5 * (qrcd(i,k) + qrcd(i,k+1))
            xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo
c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(dwnflg2(i)) then
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            if(gcmflg.and.edtx(i).ge..90) then
               cnvflg(i) = .false.
               dwnflg2(i) = .false.
            endif
            edtx(i) = min(edtx(i),edtmax)
          endif
        else
          edtx(i) = 0.
        endif
      enddo
c
c
c--- downdraft cloudwork functions
c
      do k = kmax-1, 1, -1
        do i = 1, im
          if(dwnflg2(i).and.k.lt.jmin(i)) then
            gamma = el2orc * qeso(i,k+1) / to(i,k+1)**2
            dhh=xhcd(i)
            dt= to(i,k)
            dg= gamma
            dh= heso(i,k)
            dz=-1.*(zo(i,k+1)-zo(i,k))
            xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))
     &              *(1.+fv*cp*dg*dt/hvap)
            xaa0(i)=xaa0(i)+edtx(i)*
     &      dz*g*fv*max(0.,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
c
c  calculate critical cloud work function
c
      do i = 1, im
        acrt(i) = 0.
        if(cnvflg(i)) then
          if(p(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.-p(i,ktcon(i)))
     &              /(975.-pcrit(15))
          else if(p(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k = int((850. - p(i,ktcon(i)))/50.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*
     *           (p(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
           endif
         endif
       enddo
      do i = 1, im
        acrtfct(i) = 1.
        w1 = w1s
        w2 = w2s
        w3 = w3s
        w4 = w4s
        if(slimsk(i).eq.1.) then
          w1 = w1l
          w2 = w2l
          w3 = w3l
          w4 = w4l
        endif
        if(cnvflg(i)) then
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          acrtfct(i) = max(acrtfct(i),-1.)
          acrtfct(i) = min(acrtfct(i),1.)
          acrtfct(i) = 1. - acrtfct(i)
          dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
c
c--- large scale forcing
c
      do i= 1, im
        flg(i) = cnvflg(i)
        if(cnvflg(i)) then
          f(i) = (aa1(i) - acrt(i) * acrtfct(i)) / dtconv(i)
          if(f(i).le.0.) flg(i) = .false.
        endif
        cnvflg(i) = flg(i)
        if(cnvflg(i)) then
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) flg(i) = .false.
        endif
c
c--- kernel, cloud base mass flux
c
        cnvflg(i) = flg(i)
        if(cnvflg(i)) then
          xmb(i) = -f(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
c
c  restore t0 and qo to t1 and q1 in case convection stops
c
      do k = 1, kmax
        do i = 1, im
          to(i,k) = t1(i,k)
          qo(i,k) = q1(i,k)
#ifdef ICE
          qeso(i,k) = 10. * fpvs(t1(i,k))
#else
          qeso(i,k) = 10. * fpvs0(t1(i,k))
#endif
          qeso(i,k) = eps * qeso(i,k) / (p(i,k) + epsm1 * qeso(i,k))
          qeso(i,k) = max(qeso(i,k),qmin)
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c--- feedback: simply the changes from the cloud with unit mass flux
c---           multiplied by  the mass flux necessary to keep the
c---           equilibrium with the larger-scale.
c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.k.le.ktcon(i)) then
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
            t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
            q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
            dp = 100. * psfc(i) * del(i,k)
            delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
            delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
            deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
          endif
        enddo
      enddo
      do i = 1, im
        delqbar(i) = 0.
        deltbar(i) = 0.
      enddo
      do k = kmax, 1, -1
        do i = 1, im
          dellaq(i,k) = 0.
          if(cnvflg(i).and.k.le.ktcon(i)) then
            if(evaflg.and.slimsk(i).eq.0.) then
               evef = evfact * edt(i)
            else
               evef = 0.07
            endif
c           if(evaflg.and.slimsk(i).eq.0.) evef = evfact * edt(i)
            aup = 1.
            if(k.le.kb(i)) aup = 0.
            adw = 1.
            if(k.gt.jmin(i)) adw = 0.
            rn(i) = rn(i)
     &            + (aup * pwo(i,k) + adw * edto(i) * pwdo(i,k))
     &            * xmb(i) * .001 * dt2
            qrs(i,k) = aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
            qcond(i) = evef * (qo(i,k) - qeso(i,k)) / (1. + el2orc *
     &               qeso(i,k) / to(i,k)**2)
            dp = 100. * psfc(i) * del(i,k)
            if(rn(i).gt.0..and.qcond(i).le.0.) then
              qevap = -qcond(i) * (1. - exp(-.32 * sqrt(dt2 * rn(i))))
              qevap = min(qevap, rn(i)*1000.*g/dp)
              q1(i,k) = q1(i,k) + qevap
              t1(i,k) = t1(i,k) - elocp * qevap
              rn(i) = rn(i) - .001 * qevap * dp / g
              dellat = - elocp*qevap/xmb(i)/dt2
              dellaq(i,k) =  + qevap/xmb(i)/dt2
            endif
            delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
            deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
          endif
        enddo
      enddo
c
c  precipitation rate converted to actual precip
c  in unit of m instead of kg
c
      do i = 1, im
        if(cnvflg(i)) then
c
c  in the event of upper level rain evaporation and lower level downdraf
c    moistening, rn can become negative, in this case, we back out of th
c    heating and the moistening
c
          if(rn(i).le.0.) then
            rn(i) = 0.
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kuo(i) = 1
            cldwrk(i) = aa1(i)
          endif
        endif
      enddo
      do k = 1, kmax
        do i = 1, im
          if(cnvflg(i).and.rn(i).le.0.) then
             t1(i,k) = to(i,k)
             q1(i,k) = qo(i,k)
          endif
        enddo
      enddo
      return
      end
