#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LNT2 lnt2
#define LONFP lonfp
#include <define.h>
      subroutine swr95(ipts,s0,isrc,ibnd,pl,ta,wa,oa,co2,cosz,taucl,
     1            ccly,cfac,icfc,icwp,
#ifdef ICECLOUD
     2            cwp,cip,rew,rei,fice,
#endif
     1            albuvb,albuvd,albirb,albird,paer,
     2            htrc,tupfxc,tdnflx,supfxc,sdnfxc,
     3            tupfx0,supfx0,sdnfx0,
     4            sdnfvb,sdnfvd,sdnfnb,sdnfnd,solinsw)
c     4            sdnfvb,sdnfvd,sdnfnb,sdnfnd)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c
c subprogram:    swr95      computes short-wave radiative heating
c   programmer: yu-tai hou  org: w/nmc20      date: 95-02-09
c
c abstract: this code is a modified version of m.d. chous sw
c   radiation code to fit nmc mrf and climate models.  it computes
c   sw atmospheric absorption and scattering effects due to o3,
c   h2o,co2,o2,clouds, and aerosols, etc.
c   it has 4 uv+vis bands and 3 nir bands (10 k-values each).
c
c references: chou (1986, j. clim. appl.meteor.)
c   chou (1990, j. clim.), and chou (1992, j. atms. sci.)
c
c program history log:
c   94-06-12   m.d. chou, gla.
c   95-02-09   yu-tai hou      - recode for nmc models
c
c usage:       call swr95
c
c attributes:
c   language:  fortran 77
c   machine:   cray c-90
c
c input parameters:
c   s0     : solar constant
c   isrc   : flags for selecting absorbers
c            1:aerosols, 2:o2, 3:co2, 4:h2o, 5:o3
c            =0:without it,  =1: with it.
c   pl     : model level pressure in mb
c   ta     : model layer temperature in k
c   wa     : layer specific humidity in gm/gm
c   oa     : layer ozone concentration in gm/gm
c   co2    : co2 mixing ration by volumn
c   cosz   : cosine of solar zenith angle
c   taucl  : optical depth of cloud layers
c   ccly   : layer cloud fraction
c   cfac   : fraction of clear sky view at the layer interface
c   icfc   : =0 no cloud factor to weigh clear and cloudy fluxes
c            =1 use cloud factor to weigh clear and cloudy fluxes
c   icwp   : flag indicates the method used for cloud properties
c            calculations, =0 use t-p; =1 use cwc/cic.
c   cwp    : layer cloud water path (g/m**2)
c   cip    : layer cloud ice path (g/m**2)
c   rew    : layer water cloud drop effective radius (micron)
c   rei    : layer ice cloud drop effective radius
c   fice   : fraction of cloud ice content
c   albuvb : uv+vis surf direct albedo
c   albuvd : uv+vis surf diffused albedo
c   albirb : nir surf direct albedo
c   albird : nir surf diffused albedo
c   paer   : aerosol profiles (fraction)
c
c output parameter:
c   htrc   : heating rates for cloudy sky in  k/day
c   tupfxc : upward flux at toa for cloudy sky  w/m**2
c   tdnflx : dnward flux at toa for all sky  w/m**2
c   supfxc : upward flux at sfc for cloudy sky  w/m**2
c   sdnfxc : dnward flux at sfc for cloudy sky  w/m**2
c   tupfx0 : upward flux at toa for clear sky   w/m**2
c   supfx0 : upward flux at sfc for clear sky   w/m**2
c   sdnfx0 : dnward flux at sfc for clear sky   w/m**2
c   sdnfvb : downward surface vis beam flux     w/m**2
c   sdnfnb : downward surface nir beam flux     w/m**2
c   sdnfvd : downward surface vis diff flux     w/m**2
c   sdnfnd : downward surface nir diff flux     w/m**2
c   solinsw: net solar absorption               w/m**2 
c
c note:
c   for all quantities, k=1 is the top level/layer, except
c   si and sl, for which k=1 is the surface level/layer.
c
c$$$
c
#include <paramodel.h>
#include <comio.h>
#include <rdparm.h>
c ---  input
                      d i m e n s i o n
     1  pl (imbx,lp1), ta(imbx,lp1),   wa(imbx,l),   oa(imbx,l)
     2, taucl(imbx,l), ccly(imbx,l), cfac(imbx,lp1),cosz(imax)
     3, albuvb(imax),  albuvd(imax), albirb(imax), albird(imax)
     4, paer(imbx,6),  isrc(nsrc)
#ifdef ICECLOUD
     4, fice(imax,l)
     5, cwp(imax,l),   cip(imax,l), rew(imax,l), rei(imax,l)
#endif

c ---  output
                      d i m e n s i o n
     1  tupfxc(imax), supfxc(imax), sdnfxc(imax), tdnflx(imax)
     2, tupfx0(imax), supfx0(imax), sdnfx0(imax), htrc(imbx,l)
     3, sdnfvb(imax), sdnfvd(imax), sdnfnb(imax), sdnfnd(imax)
     4, sdn0vb(imax), sdn0vd(imax), sdn0nb(imax), sdn0nd(imax)
     5, solinsw(imax,l)

c ---  internal array
                      d i m e n s i o n
     1  fnet0(imbx,lp1), fnetc(imbx,lp1), htr0 (imbx,lp1)
     2, dflx0(imbx,lp1), dflxc(imax),     dp   (imbx,l)
     3, scal (imbx,l),   swh  (imbx,lp1), so2  (imbx,lp1)
     4, wh   (imbx,l),   csm  (imax),  cf0(imax),  cf1(imax)
     5, dwsfb0(imax), dwsfd0(imax), dwsfbc(imax), dwsfdc(imax)
#ifdef ICECLOUD
      real rewi(imax,l), reii(imax,l), oh(imbx,l)
#endif
      logical  daytm(imax)
      data taucrt / 0.05 /, ifpr / 0 /
c===> ... ibnd=1:use one nir band, =2:use three nir bands
c     data ibnd / 1 /  ! define in gradini
c===> ... begin here
      if (ifpr .eq. 0) then
#ifndef RSM
#ifndef NOPRINT
        if(iope) write(6,12) (isrc(i),i=1,nsrc)
  12    format(3x,'aerosol, o2, co2, h2o, o3 =',5i3)
#endif
#endif
        ifpr = 1
      end if
c
      nday = 0
      do 20 i=1,ipts
        swh (i,1) = 0.0e0
        so2 (i,1) = 0.0e0
        tdnflx(i) = s0 * cosz(i)
        tupfxc(i) = 0.0e0
        tupfx0(i) = 0.0e0
        supfxc(i) = 0.0e0
        supfx0(i) = 0.0e0
        sdnfxc(i) = 0.0e0
        sdnfx0(i) = 0.0e0
        dflxc(i)  = 0.0e0
        cf0(i)    = cfac(i,lp1)
        cf1(i)    = 1.0e0 - cf0(i)
c===> ... csm is the effective secant of the solar xenith angle
        csm   (i) = 35.0e0/(sqrt(1224.0e0*cosz(i)*cosz(i)+1.0e0))
        daytm(i) = cosz(i) .gt. 0.0e0
        if (daytm(i)) nday = nday + 1
  20  continue
      if (nday .eq. 0) then
        do 25 i=1,ipts
          sdnfvb(i) = 0.0e0
          sdnfvd(i) = 0.0e0
          sdnfnb(i) = 0.0e0
          sdnfnd(i) = 0.0e0
  25    continue
        do 30 k=1,l
        do 30 i=1,ipts
          htrc(i,k) = 0.0e0
  30    continue
        return
      end if
c
      do 40 k=1,l
      do 40 i=1,ipts
c===> ... layer thickness and pressure scaling function for
c         water vapor absorption
        dp  (i,k) = pl(i,k+1) - pl(i,k)
        scal(i,k) = dp(i,k)
     1            * (0.5e0*(pl(i,k)+pl(i,k+1))/300.0e0)**0.8e0
c===> ... scaled absorber amounts for h2o(wh,swh), unit is g/cm**2
        wh(i,k) = 1.02e0 * wa(i,k) * scal(i,k)
c    1          * exp(0.00135e0*(ta(i,k)-240.0e0))
     1          * (1.0e0 + 0.00135e0*(ta(i,k)-240.0e0))
        swh(i,k+1) = swh(i,k) + wh(i,k)
  40  continue
c
c===> ... initialize fluxes
c
      do 80 k=1,lp1
      do 80 i=1,ipts
        fnet0(i,k) = 0.0e0
        fnetc(i,k) = 0.0e0
        dflx0(i,k) = 0.0e0
  80  continue
!
      if (icfc .eq. 1) then
        do i=1,ipts
          cfac(i,lp1) = 0.0
        end do
        do k=1,l
          do i=1,ipts
            if (cf1(i) .gt. 0.0) then
              rcf1 = 1.0 / cf1(i)
              cfac(i,k) = (cfac(i,k) - cf0(i)) * rcf1
              ccly(i,k) = ccly(i,k) * rcf1
            end if
          end do
        end do
      end if
c
c     if (icwp.ne. 1) then
c       do k=1,l
c         do i=1,ipts
c           taucl(i,k) = taucl(i,k) * ccly(i,k)
c         end do
c       end do
c     else
#ifdef ICECLOUD
      if (icwp.eq. 1) then
        do k=1,l
          do i=1,ipts
!0799       ccc = ccly(i,k) * sqrt(ccly(i,k))
            ccc = ccly(i,k)
            cwp(i,k) = cwp(i,k) * ccc
            cip(i,k) = cip(i,k) * ccc
            rewi(i,k) = 1.0 / rew(i,k)
            reii(i,k) = 1.0 / rei(i,k)
          end do
        end do
      end if
#endif
c
c===> ... compute nir fluxes
c
      if (isrc(4) .eq. 1) then
        do 100 i=1,ipts
          dwsfb0(i) = 0.0e0
          dwsfd0(i) = 0.0e0
          dwsfbc(i) = 0.0e0
          dwsfdc(i) = 0.0e0
 100    continue
c
#ifdef ICECLOUD
        call soliri(ipts,wh,ta,taucl,csm,daytm,ibnd,fice,
     1             isrc(1),paer,albirb,albird,
     1             icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#else
        call solir(ipts,wh,taucl,csm,daytm,ibnd,
     1             isrc(1),paer,albirb,albird,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#endif
c
      end if
c
c===> ... save surface nir band fluxes
c
      do 110 i=1,ipts
        sdnfnb(i) = cf0(i)*dwsfb0(i) + cf1(i)*dwsfbc(i)
        sdnfnd(i) = cf0(i)*dwsfd0(i) + cf1(i)*dwsfdc(i)
 110  continue
c
c===> ... compute uv+visible fluxes
c         scaled amounts for o3(wh), unit is (cm-amt)stp for o3.
      if (isrc(5) .eq. 1) then
c
        do 125 i=1,ipts
          dwsfb0(i) = 0.0e0
          dwsfd0(i) = 0.0e0
          dwsfbc(i) = 0.0e0
          dwsfdc(i) = 0.0e0
 125    continue
c
#ifdef ICECLOUD
!               scaled amounts for o3(wh), unit : (cm-amt)stp for o3.
        xa = 1.02 * 466.7
        do k=1,l
          do i=1,ipts
            oh(i,k) = xa * oa(i,k) * dp(i,k) + 1.0e-11
          enddo
        enddo
        call soluvi(ipts,wh,oh,ta,taucl,csm,daytm,fice,
     1             isrc(1),paer,albuvb,albuvd,
     1             icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#else
        xa = 1.02e0 * 466.7e0
        do 120 k=1,l
        do 120 i=1,ipts
          wh(i,k) = xa * oa(i,k) * dp(i,k)
 120    continue
        call soluv(ipts,wh,taucl,csm,daytm,
     1             isrc(1),paer,albuvb,albuvd,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#endif
c
      end if
c
      do 130 i=1,ipts
c===> ... save surface downward vis band fluxes
        sdnfvb(i) = cf0(i)*dwsfb0(i) + cf1(i)*dwsfbc(i)
        sdnfvd(i) = cf0(i)*dwsfd0(i) + cf1(i)*dwsfdc(i)
c===> ... compute final fluxes
        tupfxc(i) = cf0(i)*tupfx0(i) + cf1(i)*tupfxc(i)
        supfxc(i) = cf0(i)*supfx0(i) + cf1(i)*supfxc(i)
        sdnfxc(i) = cf0(i)*sdnfx0(i) + cf1(i)*sdnfxc(i)
 130  continue
      do 140 k=1,lp1
      do 140 i=1,ipts
        fnetc (i,k) = cf0(i)*fnet0(i,k) + cf1(i)*fnetc(i,k)
 140  continue
c
c===> ... compute the absorption due to oxygen,chou(1990,j.climate,209-217)
c         scaled amounts for o2(o2,so2), unit is (cm-atm)stp for o2.
      if (isrc(2) .eq. 1) then
        do 150 k=1,l
        do 150 i=1,ipts
          so2(i,k+1) = so2(i,k) + 165.22e0 * scal(i,k)
 150    continue
c===> ... to2 is the broadband transmission function for oxygen
c         0.0287 is the fraction of solar flux in the o2 bands
        do 160 k=2,lp1
        do 160 i=1,ipts
          to2 = exp(-0.00027e0 * sqrt(so2(i,k) * csm(i)))
          dflx0(i,k) = 0.0287e0 * (1.0e0 - to2)
 160    continue
      end if
c
c===> ... table look-up for the absorption due to co2
c         compute scaled amounts for co2(wc,so2).
      if (isrc(3) .eq. 1) then
        xa = co2 * 789.e0
        do 170 k=1,l
        do 170 i=1,ipts
          so2(i,k+1) = so2(i,k) + xa*scal(i,k)
 170    continue
c
        call flxco2(ipts,so2,swh,csm,daytm,dflx0)
c
      end if
c
c===> ... adjust for the effect of o2 and co2 on clear sky net fluxe
c
      if (isrc(2).eq.1 .or. isrc(3).eq.1) then
clear   do 180 k=1,lp1
c       do 180 i=1,ipts
c         fnet0(i,k) = fnet0(i,k) - dflx0(i,k)
c180    continue
c
c===> ... adjust for the effect of o2 and co2 on cloud sky net fluxe
c
        do 190 k=1,l
        do 190 i=1,ipts
          if (ccly(i,k) .gt. 0.01e0)
     1      dflxc(i) = dflxc(i) + dflx0(i,k)*cfac(i,k)*ccly(i,k)
          fnetc(i,k+1) = fnetc(i,k+1) - dflxc(i)
     1                 - dflx0(i,k+1)*cfac(i,k+1)
 190    continue
c
c===> ... adjust for other fluxes
c
        do 200 i=1,ipts
          xx = dflxc(i) + cf0(i)*dflx0(i,lp1)
          sdnfx0(i) = sdnfx0(i) - dflx0(i,lp1)
          sdnfxc(i) = sdnfxc(i) - xx
          sdnfnb(i) = sdnfnb(i) - xx
 200    continue
      end if
      if (icfc .eq. 1) then
c===> ... compute final fluxes at top and surface
        do i=1,ipts
          sdnfvb(i) = cf0(i)*sdn0vb(i) + cf1(i)*sdnfvb(i)
          sdnfvd(i) = cf0(i)*sdn0vd(i) + cf1(i)*sdnfvd(i)
          sdnfnb(i) = cf0(i)*sdn0nb(i) + cf1(i)*sdnfnb(i)
          sdnfnd(i) = cf0(i)*sdn0nd(i) + cf1(i)*sdnfnd(i)
          tupfxc(i) = cf0(i)*tupfx0(i) + cf1(i)*tupfxc(i)
          supfxc(i) = cf0(i)*supfx0(i) + cf1(i)*supfxc(i)
          sdnfxc(i) = cf0(i)*sdnfx0(i) + cf1(i)*sdnfxc(i)
        enddo
        do k=1,lp1
          do i=1,ipts
            fnetc (i,k) = cf0(i)*fnet0(i,k) + cf1(i)*fnetc(i,k)
          enddo
        enddo
      end if
c
c===> ... convert flux unit to w/m**2
c
      do 210 k=1,lp1
      do 210 i=1,ipts
clear   fnet0 (i,k) = fnet0(i,k) * tdnflx(i)
        fnetc (i,k) = fnetc(i,k) * tdnflx(i)
 210  continue
      do 220 i=1,ipts
        sdnfnb(i) = sdnfnb(i) * tdnflx(i)
        sdnfnd(i) = sdnfnd(i) * tdnflx(i)
        sdnfvb(i) = sdnfvb(i) * tdnflx(i)
        sdnfvd(i) = sdnfvd(i) * tdnflx(i)
        tupfx0(i) = tupfx0(i) * tdnflx(i)
        tupfxc(i) = tupfxc(i) * tdnflx(i)
        supfx0(i) = supfx0(i) * tdnflx(i)
        supfxc(i) = supfxc(i) * tdnflx(i)
        sdnfx0(i) = sdnfx0(i) * tdnflx(i)
        sdnfxc(i) = sdnfxc(i) * tdnflx(i)
 220  continue
c
c===> ... fac is the factor for heating rates (in k/day)
c         if use k/sec, result should be devided by 86400.
c
c     fac = 3.6*24./10.03*.98
      fac = 8.4418744e0
c
      do 230 k=1,l
      do 230 i=1,ipts
clear   htr0(i,k) = (fnet0(i,k)-fnet0(i,k+1)) * fac / dp(i,k)
        htrc(i,k) = (fnetc(i,k)-fnetc(i,k+1)) * fac / dp(i,k)
        solinsw(i,k) = (fnetc(i,k)-fnetc(i,k+1))  ! net solar absorption
 230  continue
c
      return
      end
