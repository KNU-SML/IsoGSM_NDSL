#define CAL cal
#define CP cp
#define G g
#define HFUS hfus
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
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#define SBC sbc
#define T0C t0c
#include <define.h>
      subroutine progtm(ims2,imx2,kmx,
     &                  ps,u1,v1,t1,q1,sheleg,tskin,qsurf,
     &                  smc,stc,dm,soiltyp,sigmaf,canopy,
     &                  slrad,snowmt,snowev,delt,z0rl,plantr,tg3,
     &                  gflux,f10m,u10m,v10m,t2m,q2m,zsoil,
     &                  cm, ch, rb,rhscnpy,rhsmc,aim,bim,cim,
#ifndef VIC
     &                  rcl,sl1,slk1,slimsk,inistp,lat,
#else
     &                  rcl,sl1,slk1,slimsk,inistp,lat,rho,
#endif
#ifndef RASV2
     &                  drain,evap,hflx,rnet,ep,cowave,fm,fh,wind)
#else
     &                  drain,evap,hflx,rnet,ep,cowave,fm,fh,wind,
     &                  dsfc)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
c ca is the von karman constant
#include <paramodel.h>
      parameter (charnock=.014,ca=.4)
      parameter (rd=rd_,rv=rv_,cp=cp_,g=g_,sigma=sbc_)
      parameter (eps=rd/rv,hvap=hvap_,hfus=hfus_)
      parameter (rvrdm1=rvrdm1_,t0c=t0c_,epsm1=eps-1.)
      parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
      parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)
      parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
      parameter (bb2=-.1954,cc2=.009999)
      parameter (elocp=hvap/cp,dfsnow=.31,ch2o=4.2e6,csoil=1.26e6)
      parameter (scanop=2.,cfactr=.5,zbot=-3.,tgice=271.2)
      parameter (cice=1880.*917.)
      parameter (rhoh2o=1000.,convrad=cal_*1.e4/60.)
      parameter (ctfil1=.5,ctfil2=1.-ctfil1)
      parameter (rnu=1.51e-5,arnu=.135*rnu)
      integer soiltyp
      real kt1, kt2, ktsoil
      logical flag, flagsnw
c passing array
      dimension ps(imx2),u1(imx2),v1(imx2),t1(imx2),q1(imx2)
      dimension sheleg(imx2),snowmt(imx2),snowev(imx2)
      dimension cm(imx2),ch(imx2)
      dimension tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
      dimension smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
      dimension z0rl(imx2),plantr(imx2),soiltyp(imx2),gflux(imx2)
      dimension u10m(imx2),v10m(imx2),t2m(imx2),q2m(imx2)
      dimension slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
      dimension aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
      dimension f10m(imx2),drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
      dimension evap(imx2),hflx(imx2),rnet(imx2),ep(imx2)
      dimension fm(imx2),fh(imx2),wind(imx2)
#ifdef RASV2
      dimension dsfc(imx2)
      parameter(dsfcmax=30.0,face=1.0,delx=10000.0)
#endif
c
c local array
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
      dimension rs(ILOTS),psurf(ILOTS),
     *          theta1(ILOTS)
      dimension tv1(ILOTS),tvs(ILOTS)
      dimension z1(ILOTS),thv1(ILOTS)
      dimension rho(ILOTS),qs1(ILOTS)
      dimension qss(ILOTS),snowd(ILOTS)
      dimension etpfac(ILOTS),ustar(ILOTS),
     *          tsurf(ILOTS)
      dimension q0(ILOTS),cq(ILOTS)
      dimension stsoil(ILOTS,levs_),dew(ILOTS)
      dimension edir(ILOTS),et(ILOTS,levs_),
     *          ec(ILOTS)
      dimension z0max(ILOTS),ztmax(ILOTS)
      dimension dtv(ILOTS),adtv(ILOTS)
      dimension fm10(ILOTS),fh2(ILOTS),
     *          hlinf(ILOTS)
      dimension hl1(ILOTS),pm(ILOTS),
     *          ph(ILOTS)
      dimension hl110(ILOTS),hl12(ILOTS)
      dimension rcap(ILOTS),rsmall(ILOTS)
      dimension pm10(ILOTS),ph2(ILOTS)
      dimension olinf(ILOTS),rch(ILOTS)
      dimension dft0(ILOTS),
     *          t12(ILOTS),t14(ILOTS)
      dimension delta(ILOTS),flag(ILOTS),
     *          tref(ILOTS)
      dimension twilt(ILOTS),df1(ILOTS),
     *          etp(ILOTS)
      dimension kt1(ILOTS),fx(ILOTS)
      dimension gx(ILOTS),canfac(ILOTS)
      dimension smcz(ILOTS),dmdz(ILOTS)
      dimension ddz(ILOTS),dmdz2(ILOTS)
      dimension ddz2(ILOTS),df2(ILOTS),
     *          kt2(ILOTS)
      dimension xx(ILOTS),yy(ILOTS),
     *          zz(ILOTS)
      dimension dtdz2(ILOTS),dft2(ILOTS)
      dimension dtdz1(ILOTS),dft1(ILOTS),
     *          hcpct(ILOTS)
      dimension ai(ILOTS,levs_),
     *          bi(ILOTS,levs_),
     *          ci(ILOTS,levs_)
      dimension rhstc(ILOTS,levs_)
      dimension factsnw(ILOTS),z0(ILOTS)
      dimension slwd(ILOTS),flagsnw(ILOTS)
      dimension term1(ILOTS), term2(ILOTS),
     *          partlnd(ILOTS)
      dimension restar(ILOTS), rat(ILOTS)
c-mk
      dimension stcx(ILOTS,levs_)
c
      latd = 42
      lond = 11
      delt2 = delt * 2.
      im = ims2
      km = kmx
c
c     estimate sigma ** k at 2 m
c
      sig2k = 1. - 4. * g * 2. / (cp * 280.)
c
c  initialize variables. all units are supposedly m.k.s. unless specifie
c  psurf is in pascals
c  wind is wind speed, theta1 is adiabatic surface temp from level 1
c  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
c  surface
c  convert slrad to the civilized unit from langley minute-1 k-4
c  surface roughness length is converted to m from cm
c
      xrcl = sqrt(rcl)
#ifdef RASV2
      tem=face*delx*g/(2.0*(1.0-sl1))
#endif
      do i = 1, im
        psurf(i) = 1000. * ps(i)
        slwd(i) = slrad(i) * convrad
        wind(i) = xrcl * sqrt(u1(i) * u1(i) + v1(i) * v1(i))
#ifdef RASV2
     1            +max(0.0,min(tem*dsfc(i)/psurf(i),dsfcmax))
#endif
c       wind(i) = max(wind(i),1.e-6)
        wind(i) = max(wind(i),1.)
        q0(i) = max(q1(i),1.e-9)
        tsurf(i) = tskin(i)
        theta1(i) = t1(i) / slk1
c       if(slimsk(i).eq.1..and.theta1(i).lt.tsurf(i)) then
c         wind(i)=max(wind(i),1.)
c       endif
        tv1(i) = t1(i) * (1. + rvrdm1 * q0(i))
        thv1(i) = theta1(i) * (1. + rvrdm1 * q0(i))
        tvs(i) = tsurf(i) * (1. + rvrdm1 * q0(i))
        rho(i) = (sl1 * psurf(i)) / (rd * tv1(i))
#ifdef ICE
        qs1(i) = 1000. * fpvs(t1(i))
#else
        qs1(i) = 1000. * fpvs0(t1(i))
#endif
        qs1(i) = eps * qs1(i) / (sl1 * psurf(i) + epsm1 * qs1(i))
        qs1(i) = max(qs1(i), 1.e-8)
#ifdef ICE
        qss(i) = 1000. * fpvs(tsurf(i))
#else
        qss(i) = 1000. * fpvs0(tsurf(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
        rs(i) = plantr(i)
        z0(i) = .01 * z0rl(i)
        canopy(i)= max(canopy(i),0.)
        dm(i) = 1.
        factsnw(i) = 10.
        if(slimsk(i).eq.2.) factsnw(i) = 3.
c
c  snow depth in water equivalent is converted from mm to m unit
c
        snowd(i) = sheleg(i) / 1000.
        flagsnw(i) = .false.
c
c  when snow depth is less than 1 mm, a patchy snow is assumed and
c  soil is allowed to interact with the atmosphere.
c  we should eventually move to a linear combination of soil and
c  snow under the condition of patchy snow.
c
        if(snowd(i).gt..001.or.slimsk(i).eq.2) rs(i) = 0.
        if(snowd(i).gt..001) flagsnw(i) = .true.
      enddo
c dg  if(lat.eq.latd) then
c dg    i = lond
c dg    print *, ' wind,tv1,tvs,q1,qs1,snow,slimsk=',
c dg &   wind(i),tv1(i),tvs(i),q1(i),qs1(i),snowd(i),slimsk(i)
c dg    print *, ' slrad =', slrad(i)
c dg  endif
      do i = 1, im
        if(slimsk(i).eq.0.) then
          zsoil(i,1) = 0.
        elseif(slimsk(i).eq.1.) then
          zsoil(i,1) = -.10
        else
          zsoil(i,1) = -3. / km
        endif
      enddo
 100  continue
      do k = 2, km
        do i = 1, im
          if(slimsk(i).eq.0.) then
            zsoil(i,k) = 0.
          elseif(slimsk(i).eq.1.) then
            zsoil(i,k) = zsoil(i,k-1)
     &                   + (-2. - zsoil(i,1)) / (km - 1)
          else
            zsoil(i,k) = - 3. * float(k) / float(km)
          endif
        enddo
      enddo
      do i = 1, im
        z1(i) = -rd * tv1(i) * log(sl1) / g
        drain(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          et(i,k) = 0.
          rhsmc(i,k) = 0.
          aim(i,k) = 0.
          bim(i,k) = 1.
          cim(i,k) = 0.
          stsoil(i,k) = stc(i,k)
c-mk
          stcx(i,k) = stc(i,k)
        enddo
      enddo
      do i = 1, im
        edir(i) = 0.
        ec(i) = 0.
        evap(i) = 0.
        ep(i) = 0.
        snowmt(i) = 0.
        snowev(i) = 0.
        gflux(i) = 0.
        rhscnpy(i) = 0.
        fx(i) = 0.
        etpfac(i) = 0.
        canfac(i) = 0.
      enddo
c
c  compute stability dependent exchange coefficients
c
c  this portion of the code is presently suppressed
c
      do i = 1, im
c       if(inistp.eq.1.and.slimsk(i).ne.0.) then
          ustar(i) = .1 * wind(i)
c       endif
c       if(inistp.eq.1.and.slimsk(i).eq.0.) then
        if(slimsk(i).eq.0.) then
          ustar(i) = sqrt(g * z0(i) / charnock)
        endif
      enddo
c
c  compute stability indices (rb and hlinf)
c
      do i = 1, im
        z0max(i) = min(z0(i),1. * z1(i))
        ztmax(i) = z0max(i)
        if(slimsk(i).eq.0.) then
          restar(i) = ustar(i) * z0max(i) / vis
          restar(i) = max(restar(i),.000001)
          restar(i) =  log(restar(i))
          restar(i) = min(restar(i),5.)
          restar(i) = max(restar(i),-5.)
          rat(i) = aa1 + bb1 * restar(i) + cc1 * restar(i) ** 2
          rat(i) = rat(i) / (1. + bb2 * restar(i)
     &                       + cc2 * restar(i) ** 2)
          ztmax(i) = z0max(i) * exp(-rat(i))
        endif
      enddo
      do i = 1, im
        dtv(i) = thv1(i) - tvs(i)
        adtv(i) = abs(dtv(i))
        adtv(i) = max(adtv(i),.001)
        dtv(i) = sign(1.,dtv(i)) * adtv(i)
        rb(i) = g * dtv(i) * z1(i) / (.5 * (thv1(i) + tvs(i))
     &          * wind(i) * wind(i))
        rb(i) = max(rb(i),-5000.)
        fm(i) = log((z0max(i)+z1(i)) / z0max(i))
        fh(i) = log((ztmax(i)+z1(i)) / ztmax(i))
        fm10(i) = log((z0max(i)+10.) / z0max(i))
        fh2(i) = log((ztmax(i)+2.) / ztmax(i))
        hlinf(i) = rb(i) * fm(i) * fm(i) / fh(i)
      enddo
c
c  stable case
c
      do i = 1, im
        if(dtv(i).ge.0.) then
          hl1(i) = hlinf(i)
        endif
        if(dtv(i).ge.0..and.hlinf(i).gt..25) then
          hl0inf = z0max(i) * hlinf(i) / z1(i)
          hltinf = ztmax(i) * hlinf(i) / z1(i)
          aa = sqrt(1. + 4. * alpha * hlinf(i))
          aa0 = sqrt(1. + 4. * alpha * hl0inf)
          bb = aa
          bb0 = sqrt(1. + 4. * alpha * hltinf)
          pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
          ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
          fms = fm(i) - pm(i)
          fhs = fh(i) - ph(i)
          hl1(i) = fms * fms * rb(i) / fhs
        endif
      enddo
c
c  second iteration
c
      do i = 1, im
        if(dtv(i).ge.0.) then
          hl0 = z0max(i) * hl1(i) / z1(i)
          hlt = ztmax(i) * hl1(i) / z1(i)
          aa = sqrt(1. + 4. * alpha * hl1(i))
          aa0 = sqrt(1. + 4. * alpha * hl0)
          bb = aa
          bb0 = sqrt(1. + 4. * alpha * hlt)
          pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
          ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
          hl110(i) = hl1(i) * 10. / z1(i)
          aa = sqrt(1. + 4. * alpha * hl110(i))
          pm10(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
          hl12(i) = hl1(i) * 2. / z1(i)
c         aa = sqrt(1. + 4. * alpha * hl12(i))
          bb = sqrt(1. + 4. * alpha * hl12(i))
          ph2(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
        endif
      enddo
 200  continue
c
c  unstable case
c
c
c  check for unphysical obukhov length
c
      do i = 1, im
        if(dtv(i).lt.0.) then
          olinf(i) = z1(i) / hlinf(i)
          if(abs(olinf(i)).le.50. * z0max(i)) then
            hlinf(i) = -z1(i) / (50. * z0max(i))
          endif
        endif
      enddo
c
c  get pm and ph
c
      do i = 1, im
        if(dtv(i).lt.0..and.hlinf(i).ge.-.5) then
          hl1(i) = hlinf(i)
          pm(i) = (a0 + a1 * hl1(i)) * hl1(i)
     &            / (1. + b1 * hl1(i) + b2 * hl1(i) * hl1(i))
          ph(i) = (a0p + a1p * hl1(i)) * hl1(i)
     &            / (1. + b1p * hl1(i) + b2p * hl1(i) * hl1(i))
          hl110(i) = hl1(i) * 10. / z1(i)
          pm10(i) = (a0 + a1 * hl110(i)) * hl110(i)
     &            / (1. + b1 * hl110(i) + b2 * hl110(i) * hl110(i))
          hl12(i) = hl1(i) * 2. / z1(i)
          ph2(i) = (a0p + a1p * hl12(i)) * hl12(i)
     &            / (1. + b1p * hl12(i) + b2p * hl12(i) * hl12(i))
        endif
        if(dtv(i).lt.0.and.hlinf(i).lt.-.5) then
          hl1(i) = -hlinf(i)
          pm(i) = log(hl1(i)) + 2. * hl1(i) ** (-.25) - .8776
          ph(i) = log(hl1(i)) + .5 * hl1(i) ** (-.5) + 1.386
          hl110(i) = hl1(i) * 10. / z1(i)
          pm10(i) = log(hl110(i)) + 2. * hl110(i) ** (-.25) - .8776
          hl12(i) = hl1(i) * 2. / z1(i)
          ph2(i) = log(hl12(i)) + .5 * hl12(i) ** (-.5) + 1.386
        endif
      enddo
c
c  finish the exchange coefficient computation to provide fm and fh
c
      do i = 1, im
        fm(i) = fm(i) - pm(i)
        fh(i) = fh(i) - ph(i)
        fm10(i) = fm10(i) - pm10(i)
        fh2(i) = fh2(i) - ph2(i)
        cm(i) = ca * ca / (fm(i) * fm(i))
        ch(i) = ca * ca / (fm(i) * fh(i))
        cq(i) = ch(i)
        ustar(i) = sqrt(cm(i) * wind(i) * wind(i))
      enddo
c
c  update z0 over ocean
c
      if(inistp.lt.2.and.cowave.le.0.) then
        do i = 1, im
          if(slimsk(i).eq.0.) then
            z0(i) = (charnock / g) * ustar(i) ** 2
c  new implementation of z0
c           cc = ustar(i) * z0(i) / rnu
c           pp = cc / (1. + cc)
c           ff = g * arnu / (charnock * ustar(i) ** 3)
c           z0(i) = arnu / (ustar(i) * ff ** pp)
            z0(i) = min(z0(i),.1)
            z0(i) = max(z0(i),1.e-7)
            z0rl(i) = 100. * z0(i)
          endif
        enddo
      endif
c
c  rcp = rho cp ch v
c
      do i = 1, im
        rch(i) = rho(i) * cp * ch(i) * wind(i)
      enddo
c
c  sensible and latent heat flux over open water
c
      do i = 1, im
        if(slimsk(i).eq.0.) then
          evap(i) = elocp * rch(i) * (qss(i) - q0(i))
          dm(i) = 1.
          qsurf(i) = qss(i)
        endif
      enddo
c
c  compute soil/snow/ice heat flux in preparation for surface energy
c  balance calculation
c
      do i = 1, im
        gflux(i) = 0.
        if(slimsk(i).eq.1.) then
c         dft0(i) = ktsoil(smc(i,1),soiltyp(i))
          smcz(i) = .5 * (smc(i,1) + .20)
          dft0(i) = ktsoil(smcz(i),soiltyp(i))
        elseif(slimsk(i).eq.2.) then
c  df for ice is taken from maykut and untersteiner
c  df is in si unit of w k-1 m-1
          dft0(i) = 2.2
        endif
      enddo
 300  continue
      do i = 1, im
        if(slimsk(i).ne.0.) then
c         if(snowd(i).gt..001) then
          if(flagsnw(i)) then
c
c  when snow covered, ground heat flux comes from snow
c
            tflx = min(t1(i), tsurf(i))
c           gflux(i) = -dfsnow * (tflx - stsoil(i,1))
            gflux(i) = -dfsnow * (t1(i) - stsoil(i,1))
     &                 / (factsnw(i) * max(snowd(i),.001))
          else
c           gflux(i) = dft0(i) * (stsoil(i,1) - tsurf(i))
            gflux(i) = dft0(i) * (stsoil(i,1) - t1(i))
     &                 / (-.5 * zsoil(i,1))
          endif
          gflux(i) = max(gflux(i),-200.)
          gflux(i) = min(gflux(i),+200.)
        endif
      enddo
      do i = 1, im
        flag(i) = slimsk(i).ne.0.
        partlnd(i) = 1.
        if(snowd(i).gt.0..and.snowd(i).le..001) then
          partlnd(i) = 1. - snowd(i) / .001
        endif
      enddo
c
c  compute potential evaporation for land and sea ice
c
      do i = 1, im
        if(flag(i)) then
          t12(i) = t1(i) * t1(i)
          t14(i) = t12(i) * t12(i)
c
c  rcap = fnet - sigma t**4 + gflx - rho cp ch v (t1-theta1)
c
          rcap(i) = -slwd(i) - sigma * t14(i) + gflux(i)
     &              - rch(i) * (t1(i) - theta1(i))
c
c  rsmall = 4 sigma t**3 / rch + 1
c
          rsmall(i) = 4. * sigma * t1(i) * t12(i) / rch(i) + 1.
c
c  delta = l / cp * dqs/dt
c
          delta(i) = elocp * eps * hvap * qs1(i) / (rd * t12(i))
c
c  potential evapotranspiration ( watts / m**2 ) and
c  potential evaporation
c
          term1(i) = elocp * rsmall(i) * rch(i)*(qs1(i)-q0(i))
          term2(i) = rcap(i) * delta(i)
          ep(i) = (elocp * rsmall(i) * rch(i) * (qs1(i) - q0(i))
     &              + rcap(i) * delta(i))
          etp(i) = ep(i) /
     &              (rsmall(i) * (1. + rs(i) * wind(i) * ch(i))
     &              + delta(i))
          ep(i) = ep(i) / (rsmall(i) + delta(i))
        endif
      enddo
c
c  actual evaporation over land in three parts : edir, et, and ec
c
c  direct evaporation from soil, the unit goes from m s-1 to kg m-2 s-1
c
      do i = 1, im
        flag(i) = slimsk(i).eq.1..and.ep(i).gt.0.
      enddo
      do i = 1, im
        if(flag(i)) then
          tref(i) = .75 * thsat(soiltyp(i))
          twilt(i) = twlt(soiltyp(i))
          df1(i) = funcdf(smc(i,1),soiltyp(i))
          kt1(i) = funckt(smc(i,1),soiltyp(i))
          fx(i) = -2. * df1(i) * (smc(i,1) - .23) / zsoil(i,1)
     &            - kt1(i)
          fx(i) = fx(i)
          fx(i) = min(fx(i), ep(i)/hvap)
          fx(i) = max(fx(i),0.)
c
c  sigmaf is the fraction of area covered by vegetation
c
          edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i)
        endif
      enddo
c
c  transpiration from all levels of the soil
c
      do i = 1, im
        if(flag(i)) then
          canfac(i) = (canopy(i) / scanop) ** cfactr
          etpfac(i) = sigmaf(i) * etp(i)
     &           * (1. - canfac(i)) / hvap
          gx(i) = (smc(i,1) - twilt(i)) / (tref(i) - twilt(i))
          gx(i) = max(gx(i),0.)
          gx(i) = min(gx(i),1.)
          et(i,1) = (zsoil(i,1) / zsoil(i,km)) * gx(i) * etpfac(i)
     &            * partlnd(i)
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if(flag(i)) then
            gx(i) = (smc(i,k) - twilt(i)) / (tref(i) - twilt(i))
            gx(i) = max(gx(i),0.)
            gx(i) = min(gx(i),1.)
            et(i,k) =
     &                (zsoil(i,k) - zsoil(i,k-1)) / zsoil(i,km)
     &              * gx(i) * etpfac(i) * partlnd(i)
          endif
        enddo
      enddo
 400  continue
c
c  canopy re-evaporation
c
      do i = 1, im
        if(flag(i)) then
          ec(i) = sigmaf(i) * canfac(i) * ep(i) / hvap
          ec(i) = ec(i) * partlnd(i)
c         ec(i) = min(ec(i),canopy(i)/delt2)
c-mk
c         ec(i) = min(ec(i),canopy(i)/delt)
        endif
      enddo
c
c  sum up total evaporation
c
      do i = 1, im
        if(flag(i)) then
         evap(i) = edir(i) + ec(i)
        endif
      enddo
      do k = 1, km
        do i = 1, im
          if(flag(i)) then
            evap(i) = evap(i) + et(i,k)
          endif
        enddo
      enddo
c
c  return evap unit from kg m-2 s-1 to watts m-2
c
      do i = 1, im
        if(flag(i)) then
          evap(i) = min(evap(i)*hvap,ep(i))
        endif
      enddo
c dg  if(lat.eq.latd) then
c dg    i = lond
c dg    print *, 'fx, sigmaf, edir, etpfac=', fx(i)*hvap,sigmaf(i),
c dg &          edir(i)*hvap,etpfac(i)*hvap
c dg    print *, ' et =', (et(i,k)*hvap,k=1,km)
c dg    print *, ' canfac, ec, evap', canfac(i),ec(i)*hvap,evap(i)
c dg  endif
c
c  evaporation over bare sea ice
c
      do i = 1, im
c       if(slimsk(i).eq.2.and.snowd(i).le..001) then
        if(slimsk(i).eq.2.) then
          evap(i) = partlnd(i) * ep(i)
        endif
      enddo
c
c  treat downward moisture flux situation
c  (evap was preset to zero so no update needed)
c  dew is converted from kg m-2 to m to conform to precip unit
c
      do i = 1, im
        flag(i) = slimsk(i).ne.0..and.ep(i).le.0.
        dew(i) = 0.
      enddo
      do i = 1, im
        if(flag(i)) then
c-mk
c         dew(i) = -ep(i) * delt2 / (hvap * rhoh2o)
          dew(i) = -ep(i) * delt / (hvap * rhoh2o)
          evap(i) = ep(i)
          dm(i) = 1.
        endif
      enddo
c
c  snow covered land and sea ice
c
      do i = 1, im
        flag(i) = slimsk(i).ne.0..and.snowd(i).gt.0.
      enddo
c
c  change of snow depth due to evaporation or sublimation
c
c  convert evap from kg m-2 s-1 to m s-1 to determine the reduction of s
c
      do i = 1, im
        if(flag(i)) then
c-mk
c         bfact = snowd(i) / (delt2 * ep(i) / ((hvap+hfus) * rhoh2o))
          bfact = snowd(i) / (delt * ep(i) / ((hvap+hfus) * rhoh2o))
          bfact = min(bfact,1.)
c
c  the evaporation of snow
c
          if(ep(i).le.0.) bfact = 1.
          if(snowd(i).le..001) then
            evap(i) = (snowd(i)/.001)*bfact*ep(i) + evap(i)
            snowev(i) = (snowd(i)/.001)*bfact*ep(i)
          else
            evap(i) = bfact * ep(i)
            snowev(i) = bfact * ep(i)
          endif
c-mk
c         tsurf(i) = t1(i) +
c    &          (rcap(i) - gflux(i) - dfsnow * (t1(i) - stsoil(i,1))
c    &           /(factsnw(i) * max(snowd(i),.001))
c    &           + theta1(i) - t1(i)
c    &           - bfact * ep(i)) / (rsmall(i) * rch(i)
c    &           + dfsnow / (factsnw(i)* max(snowd(i),.001)))
          tsurf(i) = t1(i)+(rcap(i)-evap(i))/
     1               (rsmall(i)*rch(i)+
     2                dfsnow/(factsnw(i)*max(snowd(i),.001)))
          snowd(i) = snowd(i) - snowev(i)*delt/(rhoh2o * (hvap+hfus))
          snowd(i) = max(snowd(i),0.)
        endif
      enddo
c
c  snow melt rate (m s-1)
c
 500  continue
      do i = 1, im
        flag(i) = slimsk(i).ne.0.
     &            .and.snowd(i).gt..0
      enddo
      do i = 1, im
        if(flag(i).and.tsurf(i).gt.t0c) then
          snowmt(i) = rch(i) * rsmall(i)
     &              * (tsurf(i) - t0c) / (rhoh2o * hfus)
c-mk
          snowmt(i) = min(snowmt(i),snowd(i) / delt)
c-mk
          snowd(i) = max(snowd(i) - snowmt(i) * delt, 0.)
          tsurf(i) = max(t0c,tsurf(i)
     &             -hfus*snowmt(i)*rhoh2o/(rch(i)*rsmall(i)))
c
c  we need to re-evaluate evaporation because of snow melt
c    the skin temperature is now bounded to 0 deg c
c
#ifdef ICE
          qss(i) = 1000. * fpvs(tsurf(i))
#else
          qss(i) = 1000. * fpvs0(tsurf(i))
#endif
          qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
          evap(i) = elocp * rch(i) * (qss(i) - q0(i))
        endif
      enddo
c
c  prepare tendency terms for the soil moisture field without precipitat
c  the unit of moisture flux needs to become m s-1 for soil moisture
c   hence the factor of rhoh2o
c
      do i = 1, im
        flag(i) = slimsk(i).eq.1.
      enddo
      do i = 1, im
        if(flag(i)) then
c-mk
c         rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt2
          rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt
          smcz(i) = max(smc(i,1), smc(i,2))
          dmdz(i) = (smc(i,1) - smc(i,2)) / (-.5 * zsoil(i,2))
          df1(i) = funcdf(smcz(i),soiltyp(i))
          kt1(i) = funckt(smcz(i),soiltyp(i))
          rhsmc(i,1) = (df1(i) * dmdz(i) + kt1(i)
     &        + (edir(i) + et(i,1))) / (zsoil(i,1) * rhoh2o)
          ddz(i) = 1. / (-.5 * zsoil(i,2))
c
c  aim, bim, and cim are the elements of the tridiagonal matrix for the
c  implicit update of the soil moisture
c
          aim(i,1) = 0.
          bim(i,1) = df1(i) * ddz(i) / (-zsoil(i,1) * rhoh2o)
          cim(i,1) = -bim(i,1)
        endif
      enddo
      do k = 2, km
        if(k.lt.km) then
          do i = 1, im
            if(flag(i)) then
              dmdz2(i) = (smc(i,k) - smc(i,k+1))
     &                   / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
              smcz(i) = max(smc(i,k), smc(i,k+1))
              df2(i) = funcdf(smcz(i),soiltyp(i))
              kt2(i) = funckt(smcz(i),soiltyp(i))
              rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)
     &             - df1(i) * dmdz(i) - kt1(i) + et(i,k))
     &                     / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
              ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
              cim(i,k) = -df2(i) * ddz2(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
            endif
          enddo
        else
          do i = 1, im
            if(flag(i)) then
              kt2(i) = funckt(smc(i,k),soiltyp(i))
              rhsmc(i,k) = (kt2(i)
     &             - df1(i) * dmdz(i) - kt1(i) + et(i,k))
     &                     / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
              drain(i) = kt2(i)
              cim(i,k) = 0.
            endif
          enddo
        endif
        do i = 1, im
          if(flag(i)) then
            aim(i,k) = -df1(i) * ddz(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
            bim(i,k) = -(aim(i,k) + cim(i,k))
            df1(i) = df2(i)
            kt1(i) = kt2(i)
            dmdz(i) = dmdz2(i)
            ddz(i) = ddz2(i)
          endif
        enddo
      enddo
 600  continue
c
c  update soil temperature and sea ice temperature
c
      do i = 1, im
        flag(i) = slimsk(i).ne.0.
      enddo
c
c  surface temperature is part of the update when snow is absent
c
      do i = 1, im
c       if(flag(i).and.snowd(i).le..001) then
        if(flag(i).and..not.flagsnw(i)) then
          yy(i) = t1(i) +
c-mk
c    &          (rcap(i)-gflux(i) + theta1(i) - t1(i)
     &          (rcap(i)-gflux(i) 
     &           - evap(i)) / (rsmall(i) * rch(i))
          zz(i) = 1. + dft0(i) / (-.5 * zsoil(i,1) * rch(i) * rsmall(i))
          xx(i) = dft0(i) * (stsoil(i,1) - yy(i)) /
     &            (.5 * zsoil(i,1) * zz(i))
        endif
c       if(flag(i).and.snowd(i).gt..001) then
        if(flag(i).and.flagsnw(i)) then
          yy(i) = stsoil(i,1)
c
c  heat flux from snow is explicit in time
c
          zz(i) = 1.
          xx(i) = dfsnow * (stsoil(i,1) - tsurf(i))
     &            / (-factsnw(i) * max(snowd(i),.001))
        endif
      enddo
c
c  compute the forcing and the implicit matrix elements for update
c
c  ch2o is the heat capacity of water and csoil is the heat capacity of
c
      do i = 1, im
        if(flag(i)) then
          smcz(i) = max(smc(i,1), smc(i,2))
          dtdz1(i) = (stsoil(i,1) - stsoil(i,2)) / (-.5 * zsoil(i,2))
          if(slimsk(i).eq.1.) then
            dft1(i) = ktsoil(smcz(i),soiltyp(i))
            hcpct(i) = smc(i,1) * ch2o + (1. - smc(i,1)) * csoil
          else
            dft1(i) = dft0(i)
            hcpct(i) = cice
          endif
          dft2(i) = dft1(i)
          ddz(i) = 1. / (-.5 * zsoil(i,2))
c
c  ai, bi, and ci are the elements of the tridiagonal matrix for the
c  implicit update of the soil temperature
c
          ai(i,1) = 0.
          bi(i,1) = dft1(i) * ddz(i) / (-zsoil(i,1) * hcpct(i))
          ci(i,1) = -bi(i,1)
          bi(i,1) = bi(i,1)
     &            + dft0(i) / (.5 * zsoil(i,1) **2 * hcpct(i) * zz(i))
c         ss = dft0(i) * (stsoil(i,1) - yy(i))
c    &         / (.5 * zsoil(i,1) * zz(i))
c         rhstc(i,1) = (dft1(i) * dtdz1(i) - ss)
          rhstc(i,1) = (dft1(i) * dtdz1(i) - xx(i))
     &                 / (zsoil(i,1) * hcpct(i))
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if(slimsk(i).eq.1.) then
            hcpct(i) = smc(i,k) * ch2o + (1. - smc(i,k)) * csoil
          elseif(slimsk(i).eq.2.) then
            hcpct(i) = cice
          endif
        enddo
        if(k.lt.km) then
          do i = 1, im
            if(flag(i)) then
              dtdz2(i) = (stsoil(i,k) - stsoil(i,k+1))
     &                   / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
              smcz(i) = max(smc(i,k), smc(i,k+1))
              if(slimsk(i).eq.1.) then
                dft2(i) = ktsoil(smcz(i),soiltyp(i))
              endif
              ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
              ci(i,k) = -dft2(i) * ddz2(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
            endif
          enddo
        else
c
c  at the bottom, climatology is assumed at 2m depth for land and
c  freezing temperature is assumed for sea ice at z(i,km)
          do i = 1, im
            if(slimsk(i).eq.1.) then
              dtdz2(i) = (stsoil(i,k) - tg3(i))
     &                   / (.5 * (zsoil(i,k-1) + zsoil(i,k)) - zbot)
              dft2(i) = ktsoil(smc(i,k),soiltyp(i))
              ci(i,k) = 0.
            endif
            if(slimsk(i).eq.2.) then
              dtdz2(i) = (stsoil(i,k) - tgice)
     &                   / (.5 * zsoil(i,k-1) - .5 * zsoil(i,k))
              dft2(i) = dft1(i)
              ci(i,k) = 0.
            endif
          enddo
        endif
        do i = 1, im
          if(flag(i)) then
            rhstc(i,k) = (dft2(i) * dtdz2(i) - dft1(i) * dtdz1(i))
     &                 / ((zsoil(i,k) - zsoil(i,k-1)) * hcpct(i))
            ai(i,k) = -dft1(i) * ddz(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k)) * hcpct(i))
            bi(i,k) = -(ai(i,k) + ci(i,k))
            dft1(i) = dft2(i)
            dtdz1(i) = dtdz2(i)
            ddz(i) = ddz2(i)
          endif
        enddo
      enddo
 700  continue
c
c  solve the tri-diagonal matrix
c
      do k = 1, km
        do i = 1, im
          if(flag(i))  then
c-mk
c           rhstc(i,k) = rhstc(i,k) * delt2
c           ai(i,k) = ai(i,k) * delt2
c           bi(i,k) = 1. + bi(i,k) * delt2
c           ci(i,k) = ci(i,k) * delt2
            rhstc(i,k) = rhstc(i,k) * delt
            ai(i,k) = ai(i,k) * delt
            bi(i,k) = 1. + bi(i,k) * delt
            ci(i,k) = ci(i,k) * delt
          endif
        enddo
      enddo
c  forward elimination
      do i = 1, im
        if(flag(i)) then
          ci(i,1) = -ci(i,1) / bi(i,1)
          rhstc(i,1) = rhstc(i,1) / bi(i,1)
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if(flag(i)) then
            cc = 1. / (bi(i,k) + ai(i,k) * ci(i,k-1))
            ci(i,k) = -ci(i,k) * cc
            rhstc(i,k) = (rhstc(i,k) - ai(i,k) * rhstc(i,k-1)) * cc
          endif
        enddo
      enddo
c  backward substituttion
      do i = 1, im
        if(flag(i)) then
          ci(i,km) = rhstc(i,km)
        endif
      enddo
      do k = km-1, 1
        do i = 1, im
          if(flag(i)) then
            ci(i,k) = ci(i,k) * ci(i,k+1) + rhstc(i,k)
          endif
        enddo
      enddo
c
c  update soil and ice temperature
c
      do k = 1, km
        do i = 1, im
          if(flag(i)) then
            stsoil(i,k) = stsoil(i,k) + ci(i,k)
          endif
        enddo
      enddo
c
c  update surface temperature for snow free surfaces
c
      do i = 1, im
c       if(slimsk(i).ne.0..and.snowd(i).le..001) then
        if(slimsk(i).ne.0..and..not.flagsnw(i)) then
c-mk
c         tsurf(i) = (yy(i) + (zz(i) - 1.) * stsoil(i,1)) / zz(i)
          tsurf(i) = t1(i)+(rcap(i)-evap(i))/
     1              (rsmall(i)*rch(i)+dft0(i)/(-.5 * zsoil(i,1)))
        endif
c       if(slimsk(i).eq.2..and.snowd(i).le..001) then
        if(slimsk(i).eq.2..and..not.flagsnw(i)) then
          tsurf(i) = min(tsurf(i),t0c)
        endif
      enddo
      do k = 1, km
        do i = 1, im
          if(slimsk(i).eq.2) then
            stc(i,k) = min(stsoil(i,k),t0c)
          endif
        enddo
      enddo
c
c  time filter for soil and skin temperature
c
#ifdef DBG
        print *,'in progtm'
        call maxmin(tsurf,im,1,1,1,'tsurf')
#endif
      if(inistp.eq.0) then
        do i = 1, im
          if(slimsk(i).ne.0.) then
c-mk
c           tskin(i) = ctfil1 * tsurf(i) + ctfil2 * tskin(i)
            tskin(i) = tsurf(i)
          endif
        enddo
#ifdef DBG
        call maxmin(tskin,im,1,1,1,'tskin')
#endif
        do k = 1, km
          do i = 1, im
            if(slimsk(i).ne.0.) then
c-mk
c             stc(i,k) = ctfil1 * stsoil(i,k) + ctfil2 * stc(i,k)
              stc(i,k) = stsoil(i,k)
            endif
          enddo
        enddo
      endif
c
c  gflux calculation
c
      do i = 1, im
        flag(i) = slimsk(i).ne.0.
c    &            .and.snowd(i).gt..001
     &            .and.flagsnw(i)
      enddo
      do i = 1, im
        if(flag(i)) then
c-mk
c         gflux(i) = -dfsnow * (tskin(i) - stc(i,1))
          gflux(i) = dfsnow * (stcx(i,1) - tskin(i))
     &               /(factsnw(i) * max(snowd(i),.001))
        endif
      enddo
      do i = 1, im
c       if(slimsk(i).ne.0..and.snowd(i).le..001) then
        if(slimsk(i).ne.0..and..not.flagsnw(i)) then
c-mk
c         gflux(i) = dft0(i) * (stc(i,1) - tskin(i))
          gflux(i) = dft0(i) * (stcx(i,1) - tskin(i))
     &               / (-.5 * zsoil(i,1))
        endif
      enddo
c
c  calculate sensible heat flux
c
      do i = 1, im
        hflx(i) = rch(i) * (tskin(i) - theta1(i))
      enddo
c
c  the rest of the output
c
      do i = 1, im
        qsurf(i) = q0(i) + evap(i) / (elocp * rch(i))
        dm(i) = 1.
c
c  convert snow depth back to mm of water equivalent
c
        sheleg(i) = snowd(i) * 1000.
      enddo
#ifdef DBG
      call maxmin(qsurf,im,1,1,1,'qsurf')
      call maxmin(evap,im,1,1,1,'evap')
      call maxmin(rch,im,1,1,1,'rch')
#endif
c     if(inistp.eq.3) then
        do i = 1, im
          f10m(i) = fm10(i) / fm(i)
          u10m(i) = f10m(i) * xrcl * u1(i)
          v10m(i) = f10m(i) * xrcl * v1(i)
           t2m(i) = tskin(i) * (1. - fh2(i) / fh(i))
     &           + theta1(i) * fh2(i) / fh(i)
           t2m(i) = t2m(i) * sig2k
c          q2m(i) = qsurf(i) * (1. - fh2(i) / fh(i))
c    &           + q0(i) * fh2(i) / fh(i)
c         t2m(i) = t1(i)
c         q2m(i) = q0(i)
          if(evap(i).ge.0.) then
c
c  in case of evaporation, use the inferred qsurf to deduce q2m
c
            q2m(i) = qsurf(i) * (1. - fh2(i) / fh(i))
     &           + q0(i) * fh2(i) / fh(i)
          else
c
c  for dew formation situation, use saturated q at tskin
c
#ifdef ICE
            qss(i) = 1000. * fpvs(tskin(i))
#else
            qss(i) = 1000. * fpvs0(tskin(i))
#endif
            qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
            q2m(i) = qss(i) * (1. - fh2(i) / fh(i))
     &           + q0(i) * fh2(i) / fh(i)
          endif
#ifdef ICE
          qss(i) = 1000. * fpvs(t2m(i))
#else
          qss(i) = 1000. * fpvs0(t2m(i))
#endif
          qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
          q2m(i) = min(q2m(i),qss(i))
        enddo
c     endif
 800  continue
      do i = 1, im
        rnet(i) = -slwd(i) - sigma * tskin(i) **4
      enddo
      return
      end
