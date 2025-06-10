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
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define NSTYPE nstype
#define NVTYPE nvtype
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#define SBC sbc
#define T0C t0c
#include <define.h>
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
      subroutine progtm99(ims2,imx2,kmx,
     &                    ps,u1,v1,t1,q1,sheleg,tskin,qsurf,
     &                    smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &                    dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &                    gflux,f10m,u10m,v10m,t2m,q2m,zsoil,
     &                    cm, ch, rb,rhscnpy,rhsmc,aim,bim,cim,
#ifndef VIC
     &                    rcl,sl1,slk1,slimsk,inistp,lat,
#else
     &                    rcl,sl1,slk1,slimsk,inistp,lat,rho,
#endif
#ifndef RASV2
#ifndef HYDRO
     &                  drain,evap,hflx,rnet,ep,fm,fh,ustar,wind)
#else
     &                  drain,evap,hflx,rnet,ep,fm,fh,ustar,wind,
     &                  ec,hydrow)
#endif
#else
     &                  drain,evap,hflx,rnet,ep,fm,fh,ustar,wind,
#ifndef HYDRO
     &                  dsfc)
#else
     &                  ec,hydrow,dsfc)
#endif
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
#ifndef HYDRO
      parameter (scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2)
#else
      parameter (scanop=2.0,cfactr=.5,zbot=-3.,tgice=271.2)
#endif
      parameter (cice=1880.*917.,topt=298.)
      parameter (rhoh2o=1000.,convrad=cal_*1.e4/60.)
      parameter (ctfil1=.5,ctfil2=1.-ctfil1)
      parameter (rnu=1.51e-5,arnu=.135*rnu)
#ifndef ICE
      parameter(nx=7501)
      dimension tbpvs0(nx)
      common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
      integer soiltyp, vegtype
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
      dimension stsoil(ILOTS,lsoil_),dew(ILOTS)
      dimension edir(ILOTS),et(ILOTS,lsoil_),
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
      dimension dft0(ILOTS),t12(ILOTS),
     *          t14(ILOTS)
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
      dimension ai(ILOTS,lsoil_),
     *          bi(ILOTS,lsoil_)
      dimension ci(ILOTS,lsoil_)
      dimension rhstc(ILOTS,lsoil_)
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
      dimension stcx(ILOTS,lsoil_)
#ifdef HYDRO
      dimension hydrow(ILOTS,lsoil_)
      dimension ecx(ILOTS)
#endif
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
c 10  ...  dwarf trees and shrubs with ground cover (tundra)
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
        wind(i) = xrcl * sqrt(u1(i) * u1(i) + v1(i) * v1(i))
#ifdef RASV2
     1            +max(0.0,min(tem*dsfc(i)/psurf(i),dsfcmax))
#endif
        wind(i) = max(wind(i),1.)
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
        xj1=min(max(c1xpvs0+c2xpvs0*t1(i),1.),float(nx))
        jx1=min(xj1,nx-1.)
        fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
        qs1(i) = 1000. * fpvs01
c       qs1(i) = 1000. * fpvs0(t1(i))
#endif
        qs1(i) = eps * qs1(i) / (sl1 * psurf(i) + epsm1 * qs1(i))
        qs1(i) = max(qs1(i), 1.e-8)
        q0(i) = min(qs1(i),q0(i))
#ifdef ICE
        qss(i) = 1000. * fpvs(tsurf(i))
#else
        xj1=min(max(c1xpvs0+c2xpvs0*tsurf(i),1.),float(nx))
        jx1=min(xj1,nx-1.)
        fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
        qss(i) = 1000. * fpvs01
c       qss(i) = 1000. * fpvs0(tsurf(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
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
#ifndef HYDRO
        canopy(i)= max(canopy(i),0.)
#endif
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
#ifdef HYDRO
          hydrow(i,k)=0.
#endif
        enddo
      enddo
      do i = 1, im
        edir(i) = 0.
        ec(i) = 0.
#ifdef HYDRO
        ecx(i)=0.
#endif
        evap(i) = 0.
        ep(i) = 0.
        snowmt(i) = 0.
        snowev(i) = 0.
        gflux(i) = 0.
        rhscnpy(i) = 0.
        fx(i) = 0.
        etpfac(i) = 0.
        canfac(i) = 0.
#ifdef HYDRO
        dew(i)=0.
#endif
      enddo
#ifdef HYDRO
c
      do i = 1, im
        if(canopy(i).lt.0.) then
          ecx(i)=canopy(i)/delt
          canopy(i)= 0.
        endif
      enddo
#endif
c
c  compute stability dependent exchange coefficients
c
c  this portion of the code is presently suppressed
c
      do i = 1, im
c         ustar(i) = .1 * wind(i)
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
c  rat taken from zeng, zhao and dickinson 1997
          rat(i) = 2.67 * restar(i) ** .25 - 2.57
          rat(i) = min(rat(i),7.)
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
        if(z0max(i).eq.0.) then
           print *,'progtm99, z0max=0. at i=',i
#ifdef MP
#ifdef RSM
           call rmpabort
#else
           call mpabort
#endif
#else
           call abort
#endif
        elseif(ztmax(i).eq.0.) then
           print *,'ztmax=0. at i=',i
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
      do i = 1, im
        if(slimsk(i).eq.0.) then
          z0(i) = (charnock / g) * ustar(i) ** 2
c  new implementation of z0
c         cc = ustar(i) * z0(i) / rnu
c         pp = cc / (1. + cc)
c         ff = g * arnu / (charnock * ustar(i) ** 3)
c         z0(i) = arnu / (ustar(i) * ff ** pp)
          z0(i) = min(z0(i),.1)
          z0(i) = max(z0(i),1.e-7)
          z0rl(i) = 100. * z0(i)
        endif
      enddo
c     endif
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
          evap(i) = elocp * rch(i) * (qss(i) - q1(i))
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
c            tflx = min(t1(i), tsurf(i))
            gflux(i) = -dfsnow * (t1(i) - stsoil(i,1))
     &                 / (factsnw(i) * max(snowd(i),.001))
          else
            gflux(i) = dft0(i) * (stsoil(i,1) - t1(i))
     &                 / (-.5 * zsoil(i,1))
          endif
c          gflux(i) = max(gflux(i),-200.)
c          gflux(i) = min(gflux(i),+200.)
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
#ifdef HYDRO
          canopy(i)=canopy(i)-ec(i)*delt
          ec(i)=ec(i)+ecx(i)
#endif
        endif
      enddo
c
c  sum up total evaporation
c
      do i = 1, im
        if(flag(i)) then
         evap(i) = edir(i) + ec(i)
#ifdef HYDRO
         ec(i)=ec(i)*hvap
#endif
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
#ifndef HYDRO
          evap(i) = min(evap(i)*hvap,ep(i))
#else
          evap(i) = evap(i)*hvap
#endif
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
#ifndef HYDRO
        dew(i) = 0.
#endif
      enddo
      do i = 1, im
        if(flag(i)) then
c-mk
          dew(i) = -ep(i) * delt / (hvap * rhoh2o)
          evap(i) = ep(i)
          dew(i) = dew(i) * partlnd(i)
          evap(i) = evap(i) * partlnd(i)
          dm(i) = 1.
#ifdef HYDRO
          ec(i)=evap(i)*sigmaf(i)/hvap
          canopy(i)=canopy(i)-ec(i)*delt
          edir(i)=evap(i)*(1.-sigmaf(i))/hvap
          ec(i) = ec(i)*hvap
          evap(i) = edir(i)*hvap + ec(i)
#endif
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
          bfact = snowd(i) / (delt * ep(i) / ((hvap+hfus) * rhoh2o))
          bfact = min(bfact,1.)
c
c  the evaporation of snow
c
          if(ep(i).le.0.) bfact = 1.
          if(snowd(i).le..001) then
            snowev(i) = bfact * ep(i) * (1. - partlnd(i))
#ifdef HYDRO
            snoex=snowev(i)*delt/(rhoh2o * (hvap+hfus))
            if(snoex.gt.snowd(i)) then
               snowev(i)=snowd(i)*(rhoh2o * (hvap+hfus))/delt
            endif
#endif
            evap(i) = evap(i) + snowev(i)
          else
            snowev(i) = bfact * ep(i)
#ifdef HYDRO
            snoex=snowev(i)*delt/(rhoh2o * (hvap+hfus))
            if(snoex.gt.snowd(i)) then
              snowev(i)=snowd(i)*(rhoh2o * (hvap+hfus))/delt
            endif
#endif
            evap(i) = snowev(i)
          endif
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
#ifndef HYDRO
#ifdef ICE
          qss(i) = 1000. * fpvs(tsurf(i))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tsurf(i),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          qss(i) = 1000. * fpvs01
c         qss(i) = 1000. * fpvs0(tsurf(i))
#endif
          qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
          evap(i) = elocp * rch(i) * (qss(i) - q0(i))
#endif
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
#ifndef HYDRO
          rhscnpy(i) = -ec(i) + sigmaf(i) * rhoh2o * dew(i) / delt
#endif
          smcz(i) = max(smc(i,1), smc(i,2))
          dmdz(i) = (smc(i,1) - smc(i,2)) / (-.5 * zsoil(i,2))
#ifndef HYDRO
          rhsmc(i,1) = (df1(i) * dmdz(i) + kt1(i)
     &        + (edir(i) + et(i,1))) / (zsoil(i,1) * rhoh2o)
          rhsmc(i,1) = rhsmc(i,1) - (1. - sigmaf(i)) * dew(i) /
     &                 ( zsoil(i,1) * delt2)
#else
          rhsmc(i,1) = -(df1(i) * dmdz(i) + kt1(i))
     &         / (-zsoil(i,1) * rhoh2o)
          hydrow(i,1)=edir(i)+et(i,1)
#endif
          ddz(i) = 1. / (-.5 * zsoil(i,2))
c
c  aim, bim, and cim are the elements of the tridiagonal matrix for the
c  implicit update of the soil moisture
c
#ifndef HYDRO
          aim(i,1) = 0.
          bim(i,1) = df1(i) * ddz(i) / (-zsoil(i,1) * rhoh2o)
          cim(i,1) = -bim(i,1)
#else
          aim(i,1) = df1(i)*ddz(i)
          bim(i,1) = -df1(i)*ddz(i)
          cim(i,1) = -kt1(i)
#endif
        endif
      enddo
      do k = 2, km
        if(k.lt.km) then
          do i = 1, im
            if(flag(i)) then
              df2(i) = funcdf(smcz(i),soiltyp(i))
c- for turning off backgound flow, set kt2(i)=0.
#ifdef NOBSFLW
              kt2(i) = 0.
#else
              kt2(i) = funckt(smcz(i),soiltyp(i))
#endif
            endif
            if(flag(i).and.stc(i,k).lt.t0c) then
              df2(i) = 0.
              kt2(i) = 0.
            endif
            if(flag(i)) then
              dmdz2(i) = (smc(i,k) - smc(i,k+1))
     &                   / (.5 * (zsoil(i,k-1) - zsoil(i,k+1)))
              smcz(i) = max(smc(i,k), smc(i,k+1))
#ifndef HYDRO
              rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)
     &             - df1(i) * dmdz(i) - kt1(i) + et(i,k))
     &                     / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
#else
c acr   do something here?  should never get here in a 2-layer model...
              rhsmc(i,k) = (df2(i) * dmdz2(i) + kt2(i)
     &             - df1(i) * dmdz(i) - kt1(i) )
     &                     / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
              hydrow(i,k)=hydrow(i,k)+et(i,k)
#endif
              ddz2(i) = 2. / (zsoil(i,k-1) - zsoil(i,k+1))
              cim(i,k) = -df2(i) * ddz2(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
            endif
          enddo
        else
          do i = 1, im
            if(flag(i)) then
#ifdef NOBSFLW
              kt2(i) = 0.
#else
c- for turning off backgound flow, set kt2(i)=0.
              kt2(i) = funckt(smc(i,k),soiltyp(i))
#endif
            endif
            if(flag(i).and.stc(i,k).lt.t0c) kt2(i) = 0.
            if(flag(i)) then
#ifndef HYDRO
              rhsmc(i,k) = (kt2(i)
     &             - df1(i) * dmdz(i) - kt1(i) + et(i,k))
     &                     / (rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
#else
c acr  adjust rhsmc?  Not used...
              rhsmc(i,k) = (df1(i) * dmdz(i) + kt1(i))
     &                     / (-rhoh2o*(zsoil(i,k) - zsoil(i,k-1)))
              hydrow(i,k)=hydrow(i,k)+kt2(i)+et(i,k)
#endif
              drain(i) = kt2(i)
              cim(i,k) = 0.
            endif
          enddo
        endif
        do i = 1, im
#ifndef HYDRO
          if(flag(i)) then
            aim(i,k) = -df1(i) * ddz(i)
     &                / ((zsoil(i,k-1) - zsoil(i,k))*rhoh2o)
            bim(i,k) = -(aim(i,k) + cim(i,k))
            df1(i) = df2(i)
            kt1(i) = kt2(i)
            dmdz(i) = dmdz2(i)
            ddz(i) = ddz2(i)
          endif
#else
          if(flag(i)) then
            aim(i,k) = -df1(i)*ddz(i)
            bim(i,k) =  df1(i)*ddz(i)
            cim(i,k) = kt1(i)
            df1(i) = df2(i)
            kt1(i) = kt2(i)
            dmdz(i) = dmdz2(i)
            ddz(i) = ddz2(i)
          endif
#endif
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
        if(flag(i).and..not.flagsnw(i)) then
          yy(i) = t1(i) +
     &          (rcap(i)-gflux(i)
     &           - evap(i)) / (rsmall(i) * rch(i))
          zz(i) = 1. + dft0(i) / (-.5 * zsoil(i,1) * rch(i) * rsmall(i))
          xx(i) = dft0(i) * (stsoil(i,1) - yy(i)) /
     &            (.5 * zsoil(i,1) * zz(i))
        endif
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
        if(slimsk(i).ne.0..and..not.flagsnw(i)) then
          tsurf(i) = t1(i)+(rcap(i)-evap(i))/
     1              (rsmall(i)*rch(i)+dft0(i)/(-.5 * zsoil(i,1)))
        endif
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
        print *,'in progtm'
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
            stc(i,k) = stsoil(i,k)
          endif
        enddo
      enddo
c
c  gflux calculation
c
      do i = 1, im
        flag(i) = slimsk(i).ne.0.
     &            .and.flagsnw(i)
      enddo
      do i = 1, im
        if(flag(i)) then
          gflux(i) = dfsnow * (stcx(i,1) - tskin(i))
     &               /(factsnw(i) * max(snowd(i),.001))
        endif
      enddo
      do i = 1, im
        if(slimsk(i).ne.0..and..not.flagsnw(i)) then
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
      do i = 1, im
        f10m(i) = fm10(i) / fm(i)
        f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * xrcl * u1(i)
        v10m(i) = f10m(i) * xrcl * v1(i)
         t2m(i) = tskin(i) * (1. - fh2(i) / fh(i))
     &         + theta1(i) * fh2(i) / fh(i)
         t2m(i) = t2m(i) * sig2k
        if(evap(i).ge.0.) then
c
c  in case of evaporation, use the inferred qsurf to deduce q2m
c
          q2m(i) = qsurf(i) * (1. - fh2(i) / fh(i))
     &         + q1(i) * fh2(i) / fh(i)
        else
c
c  for dew formation situation, use saturated q at tskin
c
#ifdef ICE
          qss(i) = 1000. * fpvs(tskin(i))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*tskin(i),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          qss(i) = 1000. * fpvs01
c         qss(i) = 1000. * fpvs0(tskin(i))
#endif
          qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
          q2m(i) = qss(i) * (1. - fh2(i) / fh(i))
     &         + q1(i) * fh2(i) / fh(i)
        endif
#ifdef ICE
        qss(i) = 1000. * fpvs(t2m(i))
#else
        xj1=min(max(c1xpvs0+c2xpvs0*t2m(i),1.),float(nx))
        jx1=min(xj1,nx-1.)
        fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
        qss(i) = 1000. * fpvs01
c       qss(i) = 1000. * fpvs0(t2m(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
        q2m(i) = min(q2m(i),qss(i))
      enddo
 800  continue
      do i = 1, im
        rnet(i) = -slwd(i) - sigma * tskin(i) **4
      enddo
      return
      end
