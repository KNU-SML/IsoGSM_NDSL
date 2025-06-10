#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine solir(ipts,wh,taucl,csm,daytm,ibnd,
     1                 kaer,paer,albb,albd,
     2                 tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3                 fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c********************************************************************
c  compute solar flux in the nir region (3 bands, 10-k per band)
c  the nir region has three water vapor bands, ten ks for each band.
c    1.   1000-4400 (/cm)         2.27-10.0 (micron)
c    2.   4400-8200               1.22-2.27
c    3.   8200-14300              0.70-1.22
c
c  input parameters:                           units
c    wh,taucl,csm,daytm,ibnd,kaer,paer,albb,albd
c  fixed input data:
c    h2o absorption coefficient (xk)           cm**2/gm
c    k-distribution function    (hk)           fraction
c
c  the following parameters must specified by users:
c    cloud single scattering albedo (sacl)     n/d
c    cloud asymmetry factor (asycl)            n/d
c  aerosols optical parameters are obtained from calling
c    subprogram aeros
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
c********************************************************************
c
#include <paramodel.h>
#include <rdparm.h>
c --- input
                      d i m e n s i o n
     1  wh(imbx,l),  taucl(imbx,l),  csm(imax),  paer(imbx,nae)
     2, albb(imax),  albd(imax)
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
c    2, dwsfxb(imax),     dwsfxd(imax),       taurs (l)
     2, dwsfxb(imax),     dwsfxd(imax)
     3, tauto (imbx,l),   ssato (imbx,l),     asyto (imbx,l)
     4, taurs (l),        ssat1 (imbx,l),     asyt1 (imbx,l)
c    4, taut1 (imbx,l),   ssat1 (imbx,l),     asyt1 (imbx,l)
     5, tauaer(imbx,l),   ssaaer(imbx,l),     asyaer(imbx,l)
     6, r0  (imbx,lp1),   t0  (imbx,lp1),     tb  (imbx,lp1)
     7, rf  (imbx,lp1),   tf  (imbx,lp1)
     8, xk  (nk0),        hk  (nk0,nrb),      sacl(nrb)
c
      data xk / 0.0010, 0.0133, 0.0422, 0.1334, 0.4217,
     1          1.3340, 5.6230, 31.620, 177.80, 1000.0 /
      data hk / .01074, .00360, .00411, .00421, .00389,
     1          .00326, .00499, .00465, .00245, .00145,
     2          .08236, .01157, .01133, .01143, .01240,
     2          .01258, .01381, .00650, .00244, .00094,
     3          .20673, .03497, .03011, .02260, .01336,
     3          .00696, .00441, .00115, .00026, .00000,
     4          .29983, .05014, .04555, .03824, .02965,
     4          .02280, .02321, .01230, .00515, .00239 /
c
      data sacl / 0.98, 0.994, 0.9995, 0.99 /, asycl / 0.843 /
      data fpmin,fpmax /1.0e-6, 0.999999/
c     data ssawv/0.00001/
c
c===> ... loop over three nir bands
c
      if (ibnd .eq. 1) then
        ibb1 = nrb
        ibb2 = nrb
      else
        ibb1 = 1
        ibb2 = nrb - 1
      end if
      do 200 ib=ibb1,ibb2
c
c===> ... get aerosols and rayleigh scattering optical properties
c
      ib1 = nvb + ib
      call aeros(ipts,ib1,kaer,paer,tauaer,ssaaer,asyaer,taurs)
c
      ssacl=sacl(ib)
c
c===> ... ik is the index for the k-distribution function (or the
c     absorption coefficient)
c
      do 100 ik=1,nk0
c
        if (hk(ik,ib) .lt. 0.00001) go to 100
c
c===> ... compute tatal optical thickness, single scattering albedo,
c         and asymmetry factor for clear sky
c
        do 30 k=1,l
        do 30 i=1,ipts
          tauwv = xk(ik)*wh(i,k)
c         taut1(i,k) = tauwv + tauaer(i,k) + taurs(k)
          tauto(i,k) = amax1(fpmin,
     1                       tauwv+tauaer(i,k)+taurs(k))
c         ssat1(i,k) = ssawv*tauwv+ssaaer(i,k)*tauaer(i,k)+taurs(k)
          ssat1(i,k) = ssaaer(i,k)*tauaer(i,k)+taurs(k)
          asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
c         tauto(i,k) = amax1(fpmin, taut1(i,k))
          ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
          asyto(i,k) = asyt1(i,k) / amax1(fpmin, ssat1(i,k))
  30    continue
c
c===> ... clear sky fluxes calculations
c
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              2,upflux,dwflux,dwsfxb,dwsfxd)
c
        do 40 k=1,lp1
        do 40 i=1,ipts
          fnet0 (i,k) = fnet0 (i,k)
     1                + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
  40    continue
        do 50 i=1,ipts
          tupfx0(i) = tupfx0(i) + upflux(i,1)  *hk(ik,ib)
          supfx0(i) = supfx0(i) + upflux(i,lp1)*hk(ik,ib)
          sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1)*hk(ik,ib)
          dwsfb0(i) = dwsfb0(i) + dwsfxb(i)*hk(ik,ib)
          dwsfd0(i) = dwsfd0(i) + dwsfxd(i)*hk(ik,ib)
  50    continue
c
c===> ... compute tatal optical thickness, single scattering albedo,
c         and asymmetry factor for cloudy sky
c
        do 60 k=1,l
        do 60 i=1,ipts
          if (taucl(i,k) .ge. 0.001e0) then
c           tauto(i,k) = taucl(i,k) + taut1(i,k)
            tauto(i,k) = taucl(i,k) + tauto(i,k)
            ssat1(i,k) = ssacl*taucl(i,k) + ssat1(i,k)
            ssato(i,k) = amin1(fpmax, ssat1(i,k)/tauto(i,k))
            asyto(i,k) = (asycl*ssacl*taucl(i,k) + asyt1(i,k))
     1                 / ssat1(i,k)
          end if
  60    continue
c
c===> ... cloudy sky fluxes calculations
c
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              2,upflux,dwflux,dwsfxb,dwsfxd)
c
        do 70 k=1,lp1
        do 70 i=1,ipts
          fnetc(i,k) = fnetc(i,k)
     1               + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
  70    continue
        do 80 i=1,ipts
          tupfxc(i) = tupfxc(i) + upflux(i,1)  *hk(ik,ib)
          supfxc(i) = supfxc(i) + upflux(i,lp1)*hk(ik,ib)
          sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1)*hk(ik,ib)
          dwsfbc(i) = dwsfbc(i) + dwsfxb(i)*hk(ik,ib)
          dwsfdc(i) = dwsfdc(i) + dwsfxd(i)*hk(ik,ib)
  80    continue
c
 100  continue
 200  continue
c
      return
      end
