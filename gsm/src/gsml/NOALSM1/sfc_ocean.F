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
      subroutine sfc_ocean(ims2,imx2,kmx,
     &                    ps,t1,q1,tskin,qsurf,
     &                    dm,cm,ch,
     &                    sl1,slk1,slimsk,inistp,lat,
     &                    evap,hflx,wind)
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
c   passing array
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
c  initialize variables. all units are supposedly m.k.s. unless specifie
c  psurf is in pascals
c  wind is wind speed, theta1 is adiabatic surface temp from level 1
c  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
c  surface
c  convert slrad to the civilized unit from langley minute-1 k-4
c  surface roughness length is converted to m from cm
c
       do i = 1, im
        if(slimsk(i).eq.0.) then
        psurf(i) = 1000. * ps(i)
        q0(i) = max(q1(i),1.e-8)
        tsurf(i) = tskin(i)
        theta1(i) = t1(i) / slk1
        tv1(i) = t1(i) * (1. + rvrdm1 * q0(i))
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
        endif
       enddo
c
c  rcp = rho cp ch v
c
      do i = 1, im
        if(slimsk(i).eq.0.) then
        rch(i) = rho(i) * cp * ch(i) * wind(i)
        endif
      enddo
c
c  sensible and latent heat flux over open water
c
      do i = 1, im
        if(slimsk(i).eq.0.) then
          evap(i) = elocp * rch(i) * (qss(i) - q1(i))
          hflx(i) = rch(i) * (tskin(i) - theta1(i))
          qsurf(i) = qss(i)
c         qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
          dm(i) = 1.
        endif
      enddo
c
      return
      end
