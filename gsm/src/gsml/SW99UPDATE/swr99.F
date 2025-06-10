#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define LEVS levs
#define LNT2 lnt2
#define LONFP lonfp
      subroutine swr99(ipts,s0,isrc,ibnd,pl,ta,wa,oa,co2,cosz,taucl,
     1            ccly,cfac,icfc,icwp,cwp,cip,rew,rei,fice,
     1            albuvb,albuvd,albirb,albird,kprf,idxc,cmix,denn,rh,
     2            htrc,tupfxc,tdnflx,supfxc,sdnfxc,
     3            tupfx0,supfx0,sdnfx0,
     4            sdnfvb,sdnfvd,sdnfnb,sdnfnd)
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
c   it has 8 uv+vis bands and 3 nir bands (10 k-values each).
c
c references: chou (1986, j. clim. appl.meteor.)
c   chou (1990, j. clim.), and chou (1992, j. atms. sci.)
c   chou and suarez (1999, nasa/tm-1999-104606,vol.15)
c
c program history log:
c   94-06-12   m.d. chou, gla.
c   95-02-09   yu-tai hou      - recode for nmc models
c   99-09-13   yu-tai hou      - updated to chou's june 99 version
c   99-10-13   yu-tai hou      - change aerosols to opac algorithm
c                                with 5 deg global data set (1998)
c
c usage:       call swr99
c
c attributes:
c   language:  fortran 90
c   machine:   cray c-90, ibm sp, sgi
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
c   kprf   : tropospheric aerosol profile type index
c   idxc,cmix
c          : aerosol component index and mixing ratio
c   denn   : aerosol number densities of 1st and 2nd layers
c   rh     : relative humidity in fraction
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
c
c note:
c   for all quantities, k=1 is the top level/layer, except
c   si and sl, for which k=1 is the surface level/layer.
c
c$$$
c
#include <paramodel.h>
#include <comio.h>
#include <rdparm99.h>
#include <aerparm.h>
c ---  input
                      d i m e n s i o n
     1  pl (imbx,lp1), ta(imbx,lp1),   wa(imbx,l),   oa(imbx,l)
     2, taucl(imbx,l), ccly(imbx,l), cfac(imbx,lp1),cosz(imax)
     3, albuvb(imax),  albuvd(imax), albirb(imax), albird(imax)
     4, rh  (imbx,l),  isrc(nsrc),   fice(imbx,l)
     5, cwp(imbx,l),   cip(imax,l), rew(imbx,l), rei(imbx,l)
     6, denn(ndn,imax),idxc(nxc,imax),kprf(imax),cmix(nxc,imax)

c ---  output
                      d i m e n s i o n
     1  tupfxc(imax), supfxc(imax), sdnfxc(imax), tdnflx(imax)
     2, tupfx0(imax), supfx0(imax), sdnfx0(imax), htrc(imbx,l)
     3, sdnfvb(imax), sdnfvd(imax), sdnfnb(imax), sdnfnd(imax)
     4, sdn0vb(imax), sdn0vd(imax), sdn0nb(imax), sdn0nd(imax)

c ---  internal array
                      d i m e n s i o n
     1  fnet0(imbx,lp1), fnetc(imbx,lp1), htr0 (imbx,l)
     2, dflx0(imbx,lp1), dflxc(imbx,lp1), dp   (imbx,l)
     3, scal (imbx,l),   swh  (imbx,lp1), so2  (imbx,lp1)
     4, wh   (imbx,l),   oh   (imbx,l),   swu  (imbx,lp1)
     5, cf0  (imax),     cf1  (imax),   snt(imax), cnt(imax)
      real rewi(imax,l), reii(imax,l)
      logical  daytm(imax)
c
#include <co2tab_sw.h>
c
      data ifpr / 0 /
c===> ... ibnd=1:use one nir band, =2:use three nir bands
c     data ibnd / 1 /  ! define in gradini
c===> ... begin here
      if (ifpr .eq. 0) then
#ifndef NOPRINT
        if(iope) write(6,12) (isrc(i),i=1,nsrc)
  12    format(3x,'aerosol, o2, co2, h2o, o3 =',5i3)
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
        dflx0(i,1)= 0.0e0
        cf0(i)    = cfac(i,lp1)
        cf1(i)    = 1.0e0 - cf0(i)
        daytm(i) = cosz(i) .gt. 0.0e0
        if (daytm(i)) then
          snt(i) = 1.0e0 / cosz(i) ! snt = secant of solar zenith angle
          nday = nday + 1
        else
          snt(i) = 1.0e3
        end if
        sdnfvb(i) = 0.0e0
        sdnfvd(i) = 0.0e0
        sdnfnb(i) = 0.0e0
        sdnfnd(i) = 0.0e0
        sdn0vb(i) = 0.0e0
        sdn0vd(i) = 0.0e0
        sdn0nb(i) = 0.0e0
        sdn0nd(i) = 0.0e0
  20  continue
      if (nday .eq. 0) then
        do 30 k=1,l
        do 30 i=1,ipts
          htrc(i,k) = 0.0e0
  30    continue
        return
      end if
c
      tfac = 0.5e0 / 300.0e0
      do 40 k=1,l
      do 40 i=1,ipts
c===> ... layer thickness and pressure scaling function for
c         water vapor absorption
        dp  (i,k) = pl(i,k+1) - pl(i,k)
        scal(i,k) = dp(i,k) * (tfac*(pl(i,k)+pl(i,k+1)))**0.8e0
c===> ... scaled absorber amounts for h2o(wh,swh), unit is g/cm**2
c       tem     = 0.00135e0*(ta(i,k)-240.0e0)
        wh(i,k) = 1.02e0 * wa(i,k) * scal(i,k)
     1          * exp(0.00135e0*(ta(i,k)-240.0e0))
c    1          * (1.0e0 + tem + 0.5e0*tem*tem) + 1.0e-11
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
      if (icwp.ne. 1) then
        do k=1,l
          do i=1,ipts
c           taucl(i,k) = taucl(i,k) * ccly(i,k)
            taucl(i,k) = taucl(i,k) * ccly(i,k)*sqrt(ccly(i,k))
          end do
        end do
      else
        do k=1,l
          do i=1,ipts
            ccc = ccly(i,k) * sqrt(ccly(i,k))
c           ccc = ccly(i,k)
            cwp(i,k) = cwp(i,k) * ccc
            cip(i,k) = cip(i,k) * ccc
            rewi(i,k) = 1.0 / rew(i,k)
            reii(i,k) = 1.0 / rei(i,k)
          end do
        end do
      end if
c
c===> ... compute nir fluxes
c
      if (isrc(4) .eq. 1) then
c
        call solir99(ipts,wh,ta,taucl,cosz,snt,daytm,ibnd,fice,
     1             isrc(1),kprf,idxc,cmix,denn,rh,albirb,albird,
     1             icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,sdn0nb,sdn0nd,sdnfnb,sdnfnd)
c
      end if
c
c===> ... compute uv+visible fluxes
c         scaled amounts for o3(wh), unit is (cm-amt)stp for o3.
      if (isrc(5) .eq. 1) then
        xa = 1.02 * 466.7
        do k=1,l
          do i=1,ipts
            oh(i,k) = xa * oa(i,k) * dp(i,k) + 1.0e-11
          enddo
        enddo
        call soluv99(ipts,wh,oh,ta,taucl,cosz,snt,daytm,fice,
     1             isrc(1),kprf,idxc,cmix,denn,rh,albuvb,albuvd,
     1             icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2             tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3             fnet0,fnetc,sdn0vb,sdn0vd,sdnfvb,sdnfvd)
c
      end if
c
c===> ... compute the absorption due to oxygen,chou(1990,j.climate,209-217)
c         scaled amounts for o2(o2,so2), unit is (cm-atm)stp for o2.
c         the constant 165.22=(1000/980)*23.14%*(22400/32)
      if (isrc(2) .eq. 1) then
        do i=1,ipts
          cnt(i) = 165.22e0 * snt(i)
        end do
        do 150 k=1,l
        do 150 i=1,ipts
          so2(i,k+1) = so2(i,k) + cnt(i) * scal(i,k)
 150    continue
c
c===> ... compute flux reduction due to oxygen, the constant 0.0633 is
c         the fraction of insolation contained in the oxygen bands.
c         to2 is the broadband transmission function for oxygen
        do k=2,lp1
        do i=1,ipts
          to2 = exp( -0.145e-3 * sqrt(so2(i,k)) )
          dflx0(i,k) = 0.0633e0 * (1.0e0 - to2)
        end do
        end do
      end if
c
c===> ... table look-up for the absorption due to co2
c         compute scaled amounts for co2(wc,so2).
c         the constant 789=(1000/980)*(44/28.97)*(22400/44)
      if (isrc(3) .eq. 1) then
        do i=1,ipts
          cnt(i)   = co2 * snt(i)
          so2(i,1) = max(so2(i,1), 1.0e-11)
        end do
        do 170 k=1,l
        do 170 i=1,ipts
          so2(i,k+1) = so2(i,k) + 789.0 *cnt(i)*scal(i,k)
 170    continue
c
c===> ... for co2 absorption in spectrum 1.220-2.270 micron
c         both water vapor and co2 absorptions are moderate
c         so2 and swh are the co2 and water vapor amounts
c         integrated from the top of the atmosphere
        u1 = -3.0
        du = 0.15
        w1 = -4.0
        dw = 0.15
        do k=2,lp1
        do i=1,ipts
            swu(i,k) = log10(so2(i,k))
            swh(i,k) = log10(swh(i,k)*snt(i))
        end do
        end do
c
c===> ... dflx0 is the updated flux reduction
c
        call flxco299(ipts,swu,u1,du,nu,swh,w1,dw,nw,cah,daytm,dflx0)
c
c===> ... for co2 absorption in spectrum 2.270-10.00 micron
c         where the co2 absorption has a large impact on the
c         heating of middle atmosphere
        u1 = 0.250e-3
        du = 0.050e-3
        w1 = -2.0
        dw = 0.05
c===> ... co2 mixing ratio is independent of space
c         swh is the logarithm of pressure
        do k=2,lp1
        do i=1,ipts
          swu(i,k) = cnt(i)
          swh(i,k) = log10(pl(i,k))
        end do
        end do
c===> ... dflx0 is the updated flux reduction
c
        call flxco299(ipts,swu,u1,du,nx,swh,w1,dw,ny,coa,daytm,dflx0)

c
      end if
c
c===> ... adjust for the effect of o2 and co2 on clear sky net fluxe
c
      if (isrc(2).eq.1 .or. isrc(3).eq.1) then
        do 180 k=1,lp1
        do 180 i=1,ipts
          fnet0(i,k) = fnet0(i,k) - dflx0(i,k)
 180    continue
c
c===> ... adjust for the effect of o2 and co2 on cloud sky net fluxe
c
        do 190 i=1,ipts
          jtop = lp1
c===> ... above clouds
          do k=1,lp1
            dflxc(i,k) = dflx0(i,k)
            if (cfac(i,k) .lt. 1.0) then
              jtop = k
              exit
            end if
          end do
c===> ... below cloud top
          if (jtop .lt. lp1) then
            do k=jtop+1,lp1
              dflxc(i,k) = dflx0(i,k) * (fnetc(i,k)/fnet0(i,k))
            end do
          end if
          do k=1,lp1
            fnetc(i,k) = fnetc(i,k) - dflxc(i,k)
          end do
 190    continue
c
c===> ... adjust for other fluxes
c
        do 200 i=1,ipts
          sdnfx0(i) = sdnfx0(i) - dflx0(i,lp1)
          sdnfxc(i) = sdnfxc(i) - dflxc(i,lp1)
          sdn0nb(i) = sdn0nb(i) - dflx0(i,lp1)
          sdnfnb(i) = sdnfnb(i) - dflxc(i,lp1)
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
c     fac = 3.6*24./10.031*.98
      fac = 8.4410328e0
c
      do 230 k=1,l
      do 230 i=1,ipts
clear   htr0(i,k) = (fnet0(i,k)-fnet0(i,k+1)) * fac / dp(i,k)
        htrc(i,k) = (fnetc(i,k)-fnetc(i,k+1)) * fac / dp(i,k)
 230  continue
c
      return
      end
