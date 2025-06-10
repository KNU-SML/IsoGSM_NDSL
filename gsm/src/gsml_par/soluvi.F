#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine soluvi(ipts,wz,oz,ta,taucl,csm,daytm,fice,
     1                 kaer,paer,albb,albd,
     1                 icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2                 tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3                 fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
#include <rdparm8.h>
c*******************************************************************
c  compute solar flux in the uv+visible region
c  the uv+visible region is grouped into 8 bands:
c    uv-c     (.175-.225);(.225-.245,.260-.280);(.245-.260);
c    uv-b     (.280-.295);(.295-.310);(.310-.320);
c    uv-a     (.320-.400);
c    par      (.400-.700)
c
c  input parameters:                            units
c    wz,oz,ta,taucl,csm,fice,kaer,paer,albb,albd
c    icwp,cwp,cip,cclv,rew,rei
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
c
c  program history log:
c   94-06-12   m.d. chou, gla.
c   95-02-09   yu-tai hou      - recode for nmc models
c   98-08-03   yu-tai hou      - updated cloud radiative properties
c            calculation. use slingo's method (jas 1989) on water
c            cloud, ebert and curry's method (jgr 1992) on ice cloud.
c   99-03-25   yu-tai hou      - updated cloud properties use the
c            most recent chou et al. data (j. clim 1998)
c   99-04-27   yu-tai hou      - updated cloud radiative property
c            calculations use linear t-adjusted method.
c   99-09-13   yu-tai hou      - updated to chou's june,1999 version
c
c********************************************************************
c
!
!     parameter (nvbb=4)
      parameter (nvbb=8)
!
c --- input
                        d i m e n s i o n
     1  oz(imax,l),   taucl(imax,l), albb(imax), albd(imax)
     2, csm(imax),  zth(imax,l),   paer(imax,nae)
     3, ta(imax,lp1),   fice(imax,l)
     3, cwp(imax,l),  cip(imax,l),   rew(imax,l), rei(imax,l)
     4, ccly(imax,l), wz(imax,l)
     5, rewi(imax,l), reii(imax,l)
      logical daytm(imax)
c --- output
                        d i m e n s i o n
     1  fnet0 (imax,lp1), dwsfb0(imax), dwsfd0(imax)
     2, fnetc (imax,lp1), dwsfbc(imax), dwsfdc(imax)
     3, tupfxc(imax),     supfxc(imax), sdnfxc(imax)
     4, tupfx0(imax),     supfx0(imax), sdnfx0(imax)
c --- temporary array
      real
     1  upflux(imax,lp1), dwflux(imax,lp1)
     2, dwsfxb(imax),     dwsfxd(imax)
     3, tauto (imax,l),   ssato (imax,l),   asyto (imax,l)
     4, taurs (l),   ssat1 (imax,l),   asyt1 (imax,l)
     5, tauaer(imax,l),   ssaaer(imax,l),   asyaer(imax,l)
     6, fffcw (imax,l),   ffft1 (imax,l),   fffto (imax,l)
     7, asycw (imax,l),   ssacw (imax,l)
c --- solar flux and absorption coefficients
     6, ss(nvbb),         ak(nvbb),         wk(nvbb)
c0499
c --- t adjusted cld property method
     8, a0w(2), a1w(2), b0w(2), b1w(2), b0i(2), b1i(2), b2i(2)
     8, a0i(2), a1i(2), c0w(2), c1w(2), c0i(2), c1i(2), c2i(2)
     9, ssaw0(2), ssai0(2), asyw0(2), asyi0(2)
     a, facw(imax,l), faci(imax,l)
!
      logical cloudy(imax)
      integer ncloud
!     logical lprnt
!
c
      data ss / 0.00057, 0.00367, 0.00083, 0.00417,
     1          0.00600, 0.00556, 0.05913, 0.39081 /
      data ak / 30.47, 187.2, 301.9, 42.83,
     1          7.090, 1.250, .0345, .0572 /
      data wk / 7*0.0e0, 0.75e-3 /
      data ssaw0 /.999998,.999998/, ssai0 /.999994,.999995/
     1     asyw0 / 0.853,  0.853 /, asyi0 / 0.7991, 0.7998/
     2,    fffrs0 / 0.1 /
      data fpmin, fpmax / 1.0e-8, 0.999999 /
c0898 - coeff for water cloud
                         d a t a
c --- t adjusted water/ice cloud coeff.
     1   a0w / 0.2807e-1,0.2798e-1 /, a1w / 0.1307e+1,0.1309e+1 /
     2,  b0w / -.1176e-6,-.1810e-6 /, c0w / 0.8276e+0,0.8272e+0 /
     3,  b1w / 0.1770e-6,0.1778e-6 /, c1w / 0.2541e-2,0.2565e-2 /
     4,  a0i / -.3011e-4,-.5975e-5 /, a1i / 0.2519e+1,0.2517e+1 /
     5,  b0i / 0.1688e-6,0.1721e-6 /, c0i / 0.7473e+0,0.7480e+0 /
     6,  b1i / 0.9936e-7,0.9177e-7 /, c1i / 0.1015e-2,0.1015e-2 /
     7,  b2i /-.1114e-10,-.1125e-10/, c2i / -.2524e-5,-.2531e-5 /
c
      save ss, ak, wk, asyw0, asyi0, ssaw0, ssai0, fffrs0, fpmin, fpmax
      save a0w,a1w,b0w,b1w,c0w,c1w,
     1     a0i,a1i,b0i,b1i,c0i,c1i,b2i,c2i
!
      do k=1,l
        do i=1,ipts
          facw(i,k) = max(0.0, min(10.0,273.15-ta(i,k)))*0.1
          faci(i,k) = max(0.0, min(30.0,263.15-ta(i,k)))/30.0
        enddo
      enddo
      cloudy(:) = .false.
!
      if (nvb .ne. nvbb) then
         print *,' nvb=',nvb,' nvbb=',nvbb,' run stopped'
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
!
      if (icwp .ne. 1) then
        do k=1,l
        do i=1,ipts
          if (taucl(i,k) .gt. 0.0) then
c0499 - t-adj prop from specified ssa and asy
            ssa1 = (1.0-fice(i,k))*(facw(i,k) *ssaw0(1)
     1                       + (1.0-facw(i,k))*ssaw0(2) )
            ssa2 =        fice(i,k) *(faci(i,k) *ssai0(1)
     1                       + (1.0-faci(i,k))*ssai0(2) )
            ssaw1 = ssa1 * taucl(i,k)
            ssaw2 = ssa2 * taucl(i,k)
            ssacw(i,k) = ssaw1 + ssaw2
            asy1 = (1.0-fice(i,k))*(facw(i,k) *asyw0(1)
     1                       + (1.0-facw(i,k))*asyw0(2) )
            asy2 =        fice(i,k) *(faci(i,k) *asyi0(1)
     1                       + (1.0-faci(i,k))*asyi0(2) )
            asyw1 = asy1 * ssaw1
            asyw2 = asy2 * ssaw2
            asycw(i,k) = asyw1 + asyw2
            fffcw(i,k) = asy1*asyw1 + asy2*asyw2
            cloudy(i) = .true.
          else
            ssacw(i,k) = 1.0
            asycw(i,k) = 0.0
            fffcw(i,k) = 0.0
          end if
        enddo
        enddo
      else
        do k=1,l
        do i=1,ipts
          if (ccly(i,k) .gt. 0.01) then
c0499 --- t-adj prop from ice/water paths
            tau1 = cwp(i,k)*(   facw(i,k) *(a0w(1)+a1w(1)*rewi(i,k))
     1                     +(1.-facw(i,k))*(a0w(2)+a1w(2)*rewi(i,k)))
            tau2 = cip(i,k)*(   faci(i,k) *(a0i(1)+a1i(1)*reii(i,k))
     1                     +(1.-faci(i,k))*(a0i(2)+a1i(2)*reii(i,k)))
            taucl(i,k) = tau1 + tau2
            ssa1 = 1.0 - (   facw(i,k) *(b0w(1)+b1w(1)*rew(i,k))
     1                 + (1.-facw(i,k))*(b0w(2)+b1w(2)*rew(i,k)) )
            ssa2 = 1.0 - (  faci(i,k) *(b0i(1)
     1                 + (b1i(1)+b2i(1)*rei(i,k))*rei(i,k))
     2                 + (1.-faci(i,k))*(b0i(2)
     3                 + (b1i(2)+b2i(2)*rei(i,k))*rei(i,k)) )
            ssaw1 = ssa1 * tau1
            ssaw2 = ssa2 * tau2
            ssacw(i,k) = ssaw1 + ssaw2
            asy1 =     facw(i,k) *(c0w(1)+c1w(1)*rew(i,k))
     1           + (1.-facw(i,k))*(c0w(2)+c1w(2)*rew(i,k))
            asy2 =     faci(i,k) *(c0i(1)
     1           + (c1i(1)+c2i(1)*rei(i,k))*rei(i,k) )
     2           + (1.-faci(i,k))*(c0i(2)
     3           + (c1i(2)+c2i(2)*rei(i,k))*rei(i,k) )
            asyw1 = asy1 * ssaw1
            asyw2 = asy2 * ssaw2
            asycw(i,k) = asyw1 + asyw2
            fffcw(i,k) = asy1*asyw1 + asy2*asyw2
            cloudy(i)  = .true.
          else
            taucl(i,k) = 0.0
            ssacw(i,k) = 1.0
            asycw(i,k) = 0.0
            fffcw(i,k) = 0.0
          end if
        enddo
        enddo
      end if
!
      ncloud = 0
      do i=1,ipts
        if (cloudy(i)) ncloud = ncloud + 1
      enddo
c
c===> ... integration over spectral bands
c
      do iv=1,nvb
c
c===> ... layer optical depth due to rayleigh scattering
c
        do k=1,l
          do i=1,ipts
            ssaaer(i,k) = 0.0
            asyaer(i,k) = 0.0
            tauaer(i,k) = 0.0
          enddo
        enddo
        call aeros8(ipts,iv,kaer,paer,tauaer,ssaaer,asyaer,taurs)
c
c===> ... compute total optical thickness, single scattering albedo,
c         and asymmetry factor for clear sky
c
        do k=1,l
        do i=1,ipts
          tauoz = ak(iv)*oz(i,k)
          tauwv = wk(iv)*wz(i,k)
          tauto(i,k) = max(fpmin, tauoz+tauwv+tauaer(i,k)+taurs(k))
          ssat1(i,k) = ssaaer(i,k)*tauaer(i,k) + taurs(k)
          asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
          ffft1(i,k) = asyaer(i,k)*asyt1(i,k) + fffrs0*taurs(k)
!
          ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
          tem        = 1.0 / max(fpmin, ssat1(i,k))
          asyto(i,k) = asyt1(i,k) * tem
          fffto(i,k) = ffft1(i,k) * tem
        enddo
        enddo
c
c===> ... clear sky fluxes calculations
c

c       call swflux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
c    1              upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              1,upflux,dwflux,dwsfxb,dwsfxd)
c
        do k=1,lp1
          do i=1,ipts
            fnet0(i,k) = fnet0(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
          enddo
        enddo
        do i=1,ipts
          tupfx0(i) = tupfx0(i) + upflux(i,1)   * ss(iv)
          supfx0(i) = supfx0(i) + upflux(i,lp1) * ss(iv)
          sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1) * ss(iv)
          dwsfb0(i) = dwsfb0(i) + dwsfxb(i)     * ss(iv)
          dwsfd0(i) = dwsfd0(i) + dwsfxd(i)     * ss(iv)
        enddo
!
        if (ncloud .gt. 0) then
c
c===> ... compute total optical thickness, single scattering albedo,
c         and asymmetry factor for cloudy sky
c
        do k=1,l
          do i=1,ipts
            if (taucl(i,k) .gt. 0.0) then
              tauto(i,k) = taucl(i,k) + tauto(i,k)
              ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
              ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
              tem        = 1.0  / max(fpmin, ssat1(i,k))
              asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
              fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
            end if
          enddo
        enddo
c
c===> ... cloudy sky fluxes calculations
c
c       call swflux(tauto,ssato,asyto,fffto,csm,zth,albb,albd,
c    1              upflux,dwflux,dwsfxb,dwsfxd, l, lp1, ipts)
        call swflux(ipts,tauto,ssato,asyto,csm,albb,albd,daytm,
     1              1,upflux,dwflux,dwsfxb,dwsfxd)
c
        do k=1,lp1
          do i=1,ipts
            fnetc(i,k) = fnetc(i,k) + (dwflux(i,k) - upflux(i,k))*ss(iv)
          enddo
        enddo
        do i=1,ipts
          tupfxc(i) = tupfxc(i) + upflux(i,1)   * ss(iv)
          supfxc(i) = supfxc(i) + upflux(i,lp1) * ss(iv)
          sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1) * ss(iv)
          dwsfbc(i) = dwsfbc(i) + dwsfxb(i)     * ss(iv)
          dwsfdc(i) = dwsfdc(i) + dwsfxd(i)     * ss(iv)
        enddo
        else
          do k=1,lp1
            do i=1,ipts
              fnetc(i,k) = fnet0(i,k)
            enddo
          enddo
          do i=1,ipts
            tupfxc(i) = tupfx0(i)
            supfxc(i) = supfx0(i)
            sdnfxc(i) = sdnfx0(i)
            dwsfbc(i) = dwsfb0(i)
            dwsfdc(i) = dwsfd0(i)
          enddo
        endif
c
      enddo           !    integration over spectral bands loop end
c
      return
      end
