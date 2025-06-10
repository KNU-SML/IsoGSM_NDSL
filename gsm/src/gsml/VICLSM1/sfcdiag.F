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
      subroutine sfcdiag(ims2,imx2,kmx,
     &                   ps,u1,v1,t1,q1,tskin,
     &                   slrad,rnet,
     &                   f10m,u10m,v10m,t2m,q2m,
     &                   rcl,slk1,slimsk,inistp,lat,
     &                   qsurf,evap,fm,fh,fm10,fh2)

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
      dimension stsoil(ILOTS,levs_),
     *          dew(ILOTS)
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
c     estimate sigma ** k at 2 m
c
      sig2k = 1. - 4. * g * 2. / (cp * 280.)
c
c  initialize variables.
c
      xrcl = sqrt(rcl)
      do i = 1, im
        psurf(i) = 1000. * ps(i)
        slwd(i) = slrad(i) * convrad
        theta1(i) = t1(i) / slk1
      enddo

c
c update surface layer properties
c
      do i = 1, im
        f10m(i) = fm10(i) / fm(i)
        f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * xrcl * u1(i)
        v10m(i) = f10m(i) * xrcl * v1(i)
         t2m(i) = tskin(i) * (1. - fh2(i) / fh(i))
     &         + theta1(i) * fh2(i) / fh(i)
         t2m(i) = t2m(i) * sig2k
c        q2m(i) = qsurf(i) * (1. - fh2(i) / fh(i))
c    &         + q1(i) * fh2(i) / fh(i)
c       t2m(i) = t1(i)
c       q2m(i) = q1(i)
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
          qss(i) = 1000. * fpvs0(tskin(i))
#endif
          qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
          q2m(i) = qss(i) * (1. - fh2(i) / fh(i))
     &         + q1(i) * fh2(i) / fh(i)
        endif
#ifdef ICE
        qss(i) = 1000. * fpvs(t2m(i))
#else
        qss(i) = 1000. * fpvs0(t2m(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
        q2m(i) = min(q2m(i),qss(i))
      enddo
c
c     compute net radiation
c
      do i = 1, im
        rnet(i) = -slwd(i) - sigma * tskin(i) **4
      enddo
c
      return
      end
