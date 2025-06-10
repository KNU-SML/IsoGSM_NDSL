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
      subroutine sfcdif(ims2,imx2,kmx,
     &                  ps,u1,v1,t1,q1,tskin,z0rl,
     &                  cm, ch, rb,
     &                  rcl,sl1,slk1,slimsk,inistp,lat,
#ifndef RASV2
     &                  fm,fh,ustar,wind,
#else
     &                  fm,fh,ustar,wind,dsfc,
#endif
#ifdef VIC
     &                  fm10,fh2,rho)
#else
     &                  fm10,fh2)
#endif
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
      latd = 23
      lond = 100
      delt2 = delt * 2.
      im = ims2
      km = kmx
c
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
        z0(i) = .01 * z0rl(i)
        z1(i) = -rd * tv1(i) * log(sl1) / g
      enddo
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
           print *,'z0max=0. at i=',i
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
c      
      return
      end
