#define CAL cal
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NSTYPE nstype
#define NVTYPE nvtype
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#define SBC sbc
#define T0C t0c
#include <define.h>
      subroutine osutm(ims2,imx2,kmx,
     &                    ps,      t1,q1,sheleg,tskin,qsurf,
     &                    smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &                    dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &                    gflux,                       zsoil,
     &                    cm, ch,    rhscnpy,rhsmc,aim,bim,cim,
     &                    rcl,sl1,slk1,slimsk,inistp,lat,
     &                  drain,evap,hflx,     ep,            wind)

#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
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
      parameter (scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2)
      parameter (cice=1880.*917.,topt=298.)
      parameter (rhoh2o=1000.,convrad=cal_*1.e4/60.)
      parameter (ctfil1=.5,ctfil2=1.-ctfil1)
      parameter (rnu=1.51e-5,arnu=.135*rnu)
      integer soiltyp, vegtype
      real kt1, kt2, ktsoil
      logical flag, flagsnw
c
c passing array
c
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
      dimension fm(imx2),fh(imx2),ustar(imx2),wind(imx2)
      dimension dlwflx(imx2)
      dimension vegtype(imx2)
#ifdef RASV2
      dimension dsfc(imx2)
      parameter(dsfcmax=30.0,face=1.0,delx=10000.0)
#endif
c
c local array
c
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
      dimension etpfac(ILOTS),tsurf(ILOTS)
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
      dimension twilt(ILOTS),df1(ILOTS)
      dimension kt1(ILOTS),fx(ILOTS)
      dimension gx(ILOTS),canfac(ILOTS)
      dimension smcz(ILOTS),dmdz(ILOTS)
      dimension ddz(ILOTS),dmdz2(ILOTS)
      dimension ddz2(ILOTS),
     *          df2(ILOTS),kt2(ILOTS)
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
      dimension term1(ILOTS),term2(ILOTS),
     *          partlnd(ILOTS)
      dimension restar(ILOTS),rat(ILOTS)
c
      dimension snet(ILOTS),smcdry(ILOTS)
      dimension rsmax(nvtype_),rgl(nvtype_)
      dimension rsmin(nvtype_),hs(nvtype_)
      dimension smdry(nstype_)
      dimension smref(nstype_),smwlt(nstype_)
c-mk
      dimension stcx(ILOTS,levs_)
c
c  the 13 vegetation types are:
c
c  1  ...  broadleave-evergreen trees (tropical forest)
c  2  ...  broadleave-deciduous trees
c  3  ...  broadleave and needle leave trees (mixed forest)
c  4  ...  needleleave-evergreen trees
c  5  ...  needleleave-deciduous trees (larch)
c  6  ...  broadleave trees with groundcover (savanna)
c  7  ...  groundcover only (perenial)
c  8  ...  broadleave shrubs with perenial groundcover
c  9  ...  broadleave shrubs with bare soil
c 10  ...  dwarf trees and shrubs with ground cover (trunda)
c 11  ...  bare soil
c 12  ...  cultivations (use parameters from type 7)
c 13  ...  glacial
c
c  the 12 vegetation type of USGS data are:
c
c  1  ... mixed farming tall grassland  (assumed to be type 12 of old veg)
c  2  ... tall/medium grassland, shrubland            (type  8 of old veg)
c  3  ... short grassland meadow and shrubland        (type  7 of old veg)
c  4  ... tundra                                      (type 10 of old veg)
c  5  ... sandy dessert                               (type 11 of old veg)
c  6  ... rocky dessert                               (type 11 of old veg)
c  7  ... tropical evergreen broadleaved forest       (type  1 of old veg)
c  8  ... evergreen forest, needleleaved forest       (type  4 of old veg)
c  9  ... medium grassland, woodland                  (type  9 of old veg)
c 10  ... deciduous forest                            (type  2 of old veg)
c 11  ... mixed deciduous and evergreen forest        (type  3 of old veg)
c 12  ... ice                                         (type 13 of old veg)
c
#ifndef USGS
      data rsmax/13*5000./
      data rsmin/150.,100.,125.,150.,100.,70.,40.,
     &           300.,400.,150.,999.,40.,999./
      data rgl/5*30.,65.,4*100.,999.,100.,999./
      data hs/41.69,54.53,51.93,47.35,47.35,54.53,36.35,
     &        3*42.00,999.,36.35,999./
#else
      data rsmax/12*5000./
      data rsmin/40.,300.,40.,150.,999.,999.,150.,
     &           150.,400.,100.,125.,999./
      data rgl/100.,100.,100.,100.,999.,999.,30.,30.,
     &         100.,30.,30.,999./
      data hs/36.35,42.00,36.35,42.00,999.,999.,41.69,
     &        47.35,42.00,54.53,51.93,999./
#endif
c
#ifdef STATSGO_SOIL
        data smref/0.236, 0.283, 0.312, 0.360, 0.360, 0.329, 0.314,
     &  0.387, 0.382, 0.338, 0.404, 0.412, 0.329, 0.0, 0.108, 0.283/
        data smwlt/0.010, 0.028, 0.047, 0.084, 0.084, 0.066, 0.067,
     &  0.120, 0.103, 0.100, 0.126, 0.138, 0.066, 0.0, 0.006, 0.028/
#else
c     data smmax/.421,.464,.468,.434,.406,.465,.404,.439,.421/
      data smdry/.07,.14,.22,.08,.18,.16,.12,.10,.07/
      data smref/.283,.387,.412,.312,.338,.382,.315,.329,.283/
      data smwlt/.029,.119,.139,.047,.010,.103,.069,.066,.029/
#endif
c
c  minvegtype and minsoiltype are defined over land only 
c  (not used over ocean)
c
      minvegtype=1
      maxvegtype=nvtype_
c
      minsoiltyp=1
      maxsoiltyp=nstype_
c
      latd = 23
      lond = 100
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
c
c  dlwflx has been given a negative sign for downward longwave
c  snet is the net shortwave flux
c
        snet(i) = -slwd(i) - dlwflx(i)
c
        q0(i) = max(q1(i),1.e-8)
        tsurf(i) = tskin(i)
        theta1(i) = t1(i) / slk1
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
        q0(i) = min(qs1(i),q0(i))
#ifdef ICE
        qss(i) = 1000. * fpvs(tsurf(i))
#else
        qss(i) = 1000. * fpvs0(tsurf(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
cx       rs(i) = plantr(i)
        rs(i) = 0.
        if(slimsk(i).eq.1) then
          if(vegtype(i).lt.minvegtype.or.vegtype(i).gt.maxvegtype) then
            print *,'illegal vegetation type'
            print *,'vegtype(i)=',vegtype(i)
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
        endif
        if(vegtype(i).gt.0) rs(i) = rsmin(vegtype(i))
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
c
c initialization is done inside sfcdrv routine
c
      do i = 1, im
        z1(i) = -rd * tv1(i) * log(sl1) / g
c       drain(i) = 0.
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
        rhscnpy(i) = 0.
        fx(i) = 0.
        etpfac(i) = 0.
        canfac(i) = 0.
      enddo
c
c  rcp = rho cp ch v
c
      do i = 1, im
        rch(i) = rho(i) * cp * ch(i) * wind(i)
      enddo
c
c  compute soil/snow/ice heat flux in preparation for surface energy
c  balance calculation
c
      do i = 1, im
        gflux(i) = 0.
        if(slimsk(i).eq.1.) then
          smcz(i) = .5 * (smc(i,1) + .20)
          if(soiltyp(i).lt.minsoiltyp.or.soiltyp(i).gt.maxsoiltyp) then
            print *,'illegal soil type'
            print *,'soiltyp(i)=',soiltyp(i)
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
c snow???
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
            gflux(i) = -dfsnow * (tflx - stsoil(i,1))
     &                 / (factsnw(i) * max(snowd(i),.001))
          else
            gflux(i) = dft0(i) * (stsoil(i,1) - tsurf(i))
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
      do i = 1, im
        snowev(i) = 0.
        if(snowd(i).gt..001) partlnd(i) = 0.
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
c soiltype???
          df1(i) = funcdf(smc(i,1),soiltyp(i))
          kt1(i) = funckt(smc(i,1),soiltyp(i))
        endif
        if(flag(i).and.stc(i,1).lt.t0c) then
          df1(i) = 0.
          kt1(i) = 0.
        endif
        if(flag(i)) then
c         tref(i) = .75 * thsat(soiltyp(i))
          tref(i) = smref(soiltyp(i))
c         twilt(i) = twlt(soiltyp(i))
          twilt(i) = smwlt(soiltyp(i))
#ifndef STATSGO_SOIL
          smcdry(i) = smdry(soiltyp(i))
#endif
c         fx(i) = -2. * df1(i) * (smc(i,1) - .23) / zsoil(i,1)
c    &            - kt1(i)
c
c  sigmaf is the fraction of area covered by vegetation
c
#ifdef NCAR_EDIR
          fx(i) = min(max((smc(i,1) - twilt(i)) 
     1        / (tref(i) - twilt(i)),0.),1.)
          edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i) * ep(i) / hvap
#else
          fx(i) = -2. * df1(i) * (smc(i,1) - smcdry(i)) / zsoil(i,1)
     &            - kt1(i)
          fx(i) = min(fx(i), ep(i)/hvap)
          fx(i) = max(fx(i),0.)
          edir(i) = fx(i) * (1. - sigmaf(i)) * partlnd(i)
#ifdef DBG
          if(lat.eq.latd) 
     1       print*,' i j edir smc sigmaf partland ',i,lat,edir(i),
     2       smc(i,1),sigmaf(i),partlnd(i)
#endif
#endif
        endif
      enddo
c
c  calculate stomatal resistance
c
      do i = 1, im
        if(flag(i)) then
c
c  resistance due to par. we use net solar flux as proxy at the present
c
          if(vegtype(i).lt.minvegtype.or.vegtype(i).gt.maxvegtype) then
            print *,'illegal vegetation type for use in rgl'
            print *,'vegtype(i)=',vegtype(i)
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
          ff = .55 * 2. * snet(i) / rgl(vegtype(i))
          rcs = (ff + rs(i)/rsmax(vegtype(i))) / (1. + ff)
          rcs = max(rcs,.0001)
          rct = 1.
          rcq = 1.
c
c  resistance due to thermal effect
c
c         rct = 1. - .0016 * (topt - theta1(i)) ** 2
c         rct = max(rct,.0001)
c
c  resistance due to humidity
c
c         rcq = 1. / (1. + hs(vegtype(i)) * (qs1(i) - q0(i)))
c         rcq = max(rcq,.0001)
c
c  compute resistance without the effect of soil moisture
c
          rs(i) = rs(i) / (rcs * rct * rcq)
        endif
      enddo
c
c  transpiration from all levels of the soil
c
      do i = 1, im
        if(flag(i)) then
          canfac(i) = (canopy(i) / scanop) ** cfactr
          etpfac(i) = sigmaf(i)
     &           * (1. - canfac(i)) / hvap
          gx(i) = (smc(i,1) - twilt(i)) / (tref(i) - twilt(i))
          gx(i) = max(gx(i),0.)
          gx(i) = min(gx(i),1.)
c
c  resistance due to soil moisture deficit
c
          rss = gx(i) * (zsoil(i,1) / zsoil(i,km))
          rss = max(rss,.0001)
          rsi = rs(i) / rss
c
c  transpiration a la monteith
c
          eth = (term1(i) + term2(i)) /
     &          (delta(i) + rsmall(i) * (1. + rsi * ch(i) * wind(i)))
          et(i,1) = etpfac(i) * eth
     &            * partlnd(i)
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if(flag(i)) then
            gx(i) = (smc(i,k) - twilt(i)) / (tref(i) - twilt(i))
            gx(i) = max(gx(i),0.)
            gx(i) = min(gx(i),1.)
c
c  resistance due to soil moisture deficit
c
          rss = gx(i) * ((zsoil(i,k) - zsoil(i,k-1)) / zsoil(i,km))
          rss = max(rss,1.e-6)
          rsi = rs(i) / rss
c
c  transpiration a la monteith
c
          eth = (term1(i) + term2(i)) /
     &          (delta(i) + rsmall(i) * (1. + rsi * ch(i) * wind(i)))
            et(i,k) = eth
     &               * etpfac(i) * partlnd(i)
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
          ec(i) = min(ec(i),canopy(i)/delt)
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
          dew(i) = dew(i) * partlnd(i)
          evap(i) = evap(i) * partlnd(i)
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
c           evap(i) = (snowd(i)/.001)*bfact*ep(i) + evap(i)
            snowev(i) = bfact * ep(i) * (1. - partlnd(i))
            evap(i) = evap(i) + snowev(i)
          else
c           evap(i) = bfact * ep(i)
            snowev(i) = bfact * ep(i)
            evap(i) = snowev(i)
          endif
c-mk
c         tsurf(i) = t1(i) +
c    &          (rcap(i) - gflux(i) - dfsnow * (t1(i) - stsoil(i,1))
c    &           /(factsnw(i) * max(snowd(i),.001))
cc   &           + theta1(i) - t1(i)
cc   &           - bfact * ep(i)) / (rsmall(i) * rch(i)
c    &           - snowev(i)) / (rsmall(i) * rch(i)
c    &           + dfsnow / (factsnw(i)* max(snowd(i),.001)))
cc        snowd(i) = snowd(i) - bfact * ep(i) * delt / (rhoh2o * hvap)
c         snowd(i) = snowd(i) - snowev(i) * delt2 / (rhoh2o * hvap)
          tsurf(i) = t1(i)+(rcap(i)-evap(i))/
     1               (rsmall(i)*rch(i)+
     2                dfsnow/(factsnw(i)*max(snowd(i),.001)))
          snowd(i) = snowd(i) - snowev(i)*delt/(rhoh2o * (hvap+hfus))
          snowd(i) = max(snowd(i),0.)
        endif
      enddo
c
c  snow melt (m s-1)
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
          snowmt(i) = min(snowmt(i),snowd(i)/delt)
          snowd(i) = snowd(i) - snowmt(i) * delt
          snowd(i) = max(snowd(i),0.)
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
          df1(i) = funcdf(smcz(i),soiltyp(i))
          kt1(i) = funckt(smcz(i),soiltyp(i))
        endif
        if(flag(i).and.stc(i,1).lt.t0c) then
          df1(i) = 0.
          kt1(i) = 0.
        endif
        if(flag(i)) then
c-mk
c         rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt2
          rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt
          smcz(i) = max(smc(i,1), smc(i,2))
          dmdz(i) = (smc(i,1) - smc(i,2)) / (-.5 * zsoil(i,2))
          rhsmc(i,1) = (df1(i) * dmdz(i) + kt1(i)
     &        + (edir(i) + et(i,1))) / (zsoil(i,1) * rhoh2o)
          rhsmc(i,1) = rhsmc(i,1) - (1. - sigmaf(i)) * dew(i) /
     &                 ( zsoil(i,1) * delt2)
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
              df2(i) = funcdf(smcz(i),soiltyp(i))
              kt2(i) = funckt(smcz(i),soiltyp(i))
            endif
            if(flag(i).and.stc(i,k).lt.t0c) then
              df2(i) = 0.
              kt2(i) = 0.
            endif
            if(flag(i)) then
              dmdz2(i) = (smc(i,k) - smc(i,k+1))
     &                   / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
              smcz(i) = max(smc(i,k), smc(i,k+1))
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
            endif
            if(flag(i).and.stc(i,k).lt.t0c) kt2(i) = 0.
            if(flag(i)) then
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
            stsoil(i,k) = min(stsoil(i,k),t0c)
          endif
        enddo
      enddo
c
c  time filter for soil and skin temperature
c
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
        print *,'in osutm'
        call maxmin(tsurf,im,1,1,1,'tsurf')
#ifdef MP
      endif
#endif
#endif
      if(inistp.eq.0) then
      do i = 1, im
        if(slimsk(i).ne.0.) then
c-mk
          tskin(i) = ctfil1 * tsurf(i) + ctfil2 * tskin(i)
          tskin(i) = tsurf(i)
        endif
      enddo
      endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
        call maxmin(tskin,im,1,1,1,'tskin')
#ifdef MP
      endif
#endif
#endif
      do k = 1, km
        do i = 1, im
          if(slimsk(i).ne.0.) then
c-mk
c           stc(i,k) = ctfil1 * stsoil(i,k) + ctfil2 * stc(i,k)
            stc(i,k) = stsoil(i,k)
          endif
        enddo
      enddo
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
        qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
        dm(i) = 1.
c
c  convert snow depth back to mm of water equivalent
c
        sheleg(i) = snowd(i) * 1000.
      enddo
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      call maxmin(qsurf,im,1,1,1,'qsurf')
      call maxmin(evap,im,1,1,1,'evap')
      call maxmin(rch,im,1,1,1,'rch')
#ifdef MP
      endif
#endif
#endif
c
      return
      end
