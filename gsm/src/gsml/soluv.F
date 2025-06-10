#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine soluv(ipts,oz,taucl,csm,daytm,
     1                 kaer,paer,albb,albd,
     2                 tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3                 fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c*******************************************************************
c  compute solar flux in the uv+visible region
c  the uv+visible region is grouped into 4 bands:
c  (.225-.285);(.175-.225,.285-.300);(.300-.325);(.325-.690)
c
c  input parameters:                            units
c    oz,taucl,csm,daytm,kaer,paer,albb,albd
c
c  output parameters:
c    fnet0  : clear sky net flux
c    fnetc  : cloudy sky net flux
c    tupfxc : cloudy sky upward flux at toa
c    supfxc : cloudy sky upward flux at sfc
c    sdnfxc : cloudy sky downward flux at sfc
c    tupfx0 : clear sky upward flux at toa
c    supfx0 : clear sky upward flux at sfc
c    sdnfx0 : clear sky downward flux at sfc
c    dwsfb0 : clear sky sfc down dir. flux
c    dwsfd0 : clear sky sfc down dif. flux
c    dwsfbc : cloudy sky sfc down dir. flux
c    dwsfdc : cloudy sky sfc down dif. flux
c
c  fixed input data:
c    fraction of solar flux contained
c       in the 8 bands (ss)                     fraction
c    rayleigh optical thickness (tauray)        /mb
c    ozone absorption coefficient (ak)          /(cm-atm)stp
c
c  the following parameters must be specified by users:
c    cloud asymmetry factor (asycl)             n/d
c  aerosol parameters are from subprogram aeros:
c********************************************************************
c
#include <paramodel.h>
#include <rdparm.h>
c --- input
                        d i m e n s i o n
     1  oz(imbx,l), taucl(imbx,l), albb(imax), albd(imax)
     2, csm(imax),  paer(imbx,nae)
      logical daytm(imax)
c --- output
                        d i m e n s i o n
     1  fnet0 (imbx,lp1), dwsfb0(imax), dwsfd0(imax)
     2, fnetc (imbx,lp1), dwsfbc(imax), dwsfdc(imax)
     3, tupfxc(imax),     supfxc(imax), sdnfxc(imax)
     4, tupfx0(imax),     supfx0(imax), sdnfx0(imax)
c --- temporary array
                        d i m e n s i o n
     1  upflux(imbx,lp1), dwflux(imbx,lp1)
c    2, dwsfxb(imax),     dwsfxd(imax),     taurs (l)
     2, dwsfxb(imax),     dwsfxd(imax)
     3, tauto (imbx,l),   ssato (imbx,l),   asyto (imbx,l)
     4, taurs (l),        ssat1 (imbx,l),   asyt1 (imbx,l)
c    4, taut1 (imbx,l),   ssat1 (imbx,l),   asyt1 (imbx,l)
     5, tauaer(imbx,l),   ssaaer(imbx,l),   asyaer(imbx,l)
     6, r0  (imbx,lp1),   t0  (imbx,lp1),   tb  (imbx,lp1)
     7, rf  (imbx,lp1),   tf  (imbx,lp1)
c --- solar flux and absorption coefficients
     6, ss(nvb),          ak(nvb)
c
      data ss / 0.00530, 0.00505, 0.01109, 0.44498 /
      data ak / 0.1805e+3, 0.267e+2, 0.199e+1, 0.050 /
      data asycl / 0.843 /
      data fpmin, fpmax / 1.0e-8, 0.9999999 /
c
c===> ... integration over spectral bands
c
      do 100 iv=1,nvb
c
c===> ... get aerosols and rayleigh scattering optical properties
c
        call aeros(ipts,iv,kaer,paer,tauaer,ssaaer,asyaer,taurs)
c
c===> ... compute total optical thickness, single scattering albedo,
c         and asymmetry factor for clear sky
c
        do 30 k=1,l
        do 30 i=1,ipts
c         taut1(i,k) = ak(iv)*oz(i,k) + tauaer(i,k) + taurs(k)
          tauto(i,k) = amax1(fpmin,
     1                       ak(iv)*oz(i,k)+tauaer(i,k)+taurs(k))
          ssat1(i,k) = ssaaer(i,k)*tauaer(i,k) + taurs(k)
          asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
c         tauto(i,k) = amax1(fpmin, taut1(i,k))
          ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
          asyto(i,k) = asyt1(i,k) / amax1(fpmin, ssat1(i,k))
  30    continue
c
c===> ... clear sky fluxes calculations
c
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              1,upflux,dwflux,dwsfxb,dwsfxd)
c
        do 40 k=1,lp1
        do 40 i=1,ipts
          fnet0(i,k) = fnet0(i,k)
     1               + (dwflux(i,k) - upflux(i,k))*ss(iv)
  40    continue
        do 50 i=1,ipts
          tupfx0(i) = tupfx0(i) + upflux(i,1)  *ss(iv)
          supfx0(i) = supfx0(i) + upflux(i,lp1)*ss(iv)
          sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1)*ss(iv)
          dwsfb0(i) = dwsfb0(i) + dwsfxb(i)*ss(iv)
          dwsfd0(i) = dwsfd0(i) + dwsfxd(i)*ss(iv)
  50    continue
c
c===> ... compute total optical thickness, single scattering albedo,
c         and asymmetry factor for cloudy sky
c
        do 60 k=1,l
        do 60 i=1,ipts
          if (taucl(i,k) .gt. 0.0e0) then
c           tauto(i,k) = taucl(i,k) + taut1(i,k)
            tauto(i,k) = taucl(i,k) + tauto(i,k)
            ssat1(i,k) = taucl(i,k) + ssat1(i,k)
            ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
            asyto(i,k) = (asycl*taucl(i,k)+asyt1(i,k))
     1                 / ssat1(i,k)
          end if
  60    continue
c
c===> ... cloudy sky fluxes calculations
c
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              1,upflux,dwflux,dwsfxb,dwsfxd)
c
        do 70 k=1,lp1
        do 70 i=1,ipts
          fnetc(i,k) = fnetc(i,k)
     1               + (dwflux(i,k) - upflux(i,k))*ss(iv)
  70    continue
        do 80 i=1,ipts
          tupfxc(i) = tupfxc(i) + upflux(i,1)  *ss(iv)
          supfxc(i) = supfxc(i) + upflux(i,lp1)*ss(iv)
          sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1)*ss(iv)
          dwsfbc(i) = dwsfbc(i) + dwsfxb(i)*ss(iv)
          dwsfdc(i) = dwsfdc(i) + dwsfxd(i)*ss(iv)
  80    continue
c
 100  continue
c
      return
      end
