#define CAL cal
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
#define PI pi
#define SBC sbc
#include <define.h>
      subroutine dcyc2(ims2,imx2,kmx,
     &                 solhr,slag,sinlab,coslab,sdec,cdec,
     &                 xlon,czmn,
     &                 sfcdlw,sfcnsw,
#ifdef VIC
     &                 sfcdsw,dswsfc,
#endif
     &                 tf,tsea,tsflw,swh,hlw,
     &                 dlwsfc,ulwsfc,slrad,tau)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    dcyc2       alters radiation for approx diurnal cycle
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-06
c           kenneth campana
c           yifeng cui  mpi       sdsc             02-09-05
c
c abstract: a diurnal cycle approximation is applied to previously
c   computed radiative fluxes and heating rates. first,the current
c   local-time value (for this particular model time step) of the
c   cosine solar zenith angle (cosz) is computed for all gaussian grid
c   points. shortwave (sw) heating rates which were computed with
c   latitudinal mean cosz in the separate radiation calculation
c   are weighted by the ratio of actual to mean cosz (see mrf model
c   documentation,1988,chapter 3,'radiative processes',authored by
c   k. campana,..). surface sw fluxes are also cosz weighted. surface
c   longwave (lw) flux from the atmosphere is altered each timestep
c   to account for diurnal changes of model temperature in the lower
c   atmosphere. lw heating rates from the separate radiation
c   computation are untouched.
c
c program history log:
c   88-05-06  joseph sela
c
c usage:    call dcyc2(lat,ssdec,solhr,colrad,czmn,sfcdlw,sfcnsw,
c                      tf,tov,slrad,swh,hlw,tau)
c   input argument list:
c     lat      - row number of gaussian latitude(n.h.).
c     ssdec    - sine of the solar declination for todays date-
c                         part of the output from radiation codes.
c     solhr    - time in hours after 00 hr greenwich.
c     colrad   - co-latitudes of gaussian grid in radians(n.h.).
c     czmn     - mean cosine solar zenith angle for all gaussian lats-
c                         part of the output from radiation codes.
c     sfcdlw   - downward lw flux at earth sfc(from radiation code)
c                in cal cm-2 min-1.
c     sfcnsw   - net sw flux at earth sfc(from radiation code using
c                czmn) in cal cm-2 min-1.
c     tf       - current value of model temperatures(w/o basic state)
c                in deg k.
c     tov      - basic state temperature for all model layers(deg k).
c     swh      - model lyr sw heating rates(from radiation code,using
c                czmn) - deg/sec.
c     hlw      - model lyr lw heating rates(from radiation code),
c                in deg/sec.
c
c   output argument list:
c     slrad    - surface net radiative flux (except lw upward flux
c                from sfc,which is added in progtn) -
c                units are cal cm-2 min-1 .
c     tau      - layer values of temperature tendency after adding
c                altered radiative heating rates-units are deg/sec.
c
c
c attributes:
c   language: fortran 200.
c   machine:  cyber 205.
c
c$$$
#include <paramodel.h>
      parameter (hsigma=sbc_,cnwatt=-cal_*1.e4/60.)
      dimension xlon(imx2),czmn(imx2),sfcdlw(imx2),sfcnsw(imx2),
     &          tf(imx2),tsea(imx2),tsflw(imx2),sinlab(imx2),
     &          swh(imx2,kmx),hlw(imx2,kmx),coslab(imx2),
     &          dlwsfc(imx2),ulwsfc(imx2),slrad(imx2),tau(imx2,kmx)
#ifdef VIC
      dimension dswsfc(imx2),sfcdsw(imx2)
#endif
c local dimension
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
#ifdef RSMVECTORIZE
      dimension xmu(igrd12p_*jgrd12p_)
#else
      dimension xmu(ILOTS)
#endif
c hmhj
      lon2=ims2
      levs=kmx
c-----------------------------------------------------------------------
c  compute cosine of solar zenith angle for both hemispheres.
      cns=pi_*(solhr-12.)/12.+slag
      do i=1,lon2
        ss=sinlab(i)*sdec
        cc=coslab(i)*cdec
        ch=cc*cos(xlon(i)+cns)
        xmu(i)=ch+ss
c       xmu(i)=(sinlab(i)*sdec)
c    1        +(coslab(i)*cdec)*cos(xlon(i)+cns)
      enddo
      do i=1,lon2
c  normalize by average value over radiation period for daytime.
        if(xmu(i).gt.0.01.and.czmn(i).gt.0.01) then
          xmu(i)=xmu(i)/czmn(i)
        else
          xmu(i)=0.
        endif
c  adjust longwave flux at surface to account for t changes in layer 1.
        sdlw=sfcdlw(i)*(tf(i)/tsflw(i))**4
#ifdef VIC
        dswsfc(i)=sfcdsw(i)*xmu(i)
#endif
c  return net surface radiative flux.
        slrad(i)=sfcnsw(i)*xmu(i)+sdlw
c  return downward and upward longwave flux at ground, respectively.
        dlwsfc(i)=sdlw*cnwatt
        ulwsfc(i)=hsigma*tsea(i)**4
      enddo
c  add radiative heating to temperature tendency
      do k=1,levs
        do i=1,lon2
          tau(i,k)=tau(i,k)+swh(i,k)*xmu(i)+hlw(i,k)
        enddo
      enddo
      return
      end
