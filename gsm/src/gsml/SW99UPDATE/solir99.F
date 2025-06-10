#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine solir99(ipts,wh,ta,taucl,zth,csm,daytm,ibnd,fice,
     1                 kaer,kprf,idxc,cmix,denn,rh,albb,albd,
     1                 icwp,cwp,cip,ccly,rew,rei,rewi,reii,
     2                 tupfxc,supfxc,sdnfxc,tupfx0,supfx0,sdnfx0,
     3                 fnet0,fnetc,dwsfb0,dwsfd0,dwsfbc,dwsfdc)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
#include <rdparm99.h>
#include <aerparm.h>
c********************************************************************
c  compute solar flux in the nir region (3 bands, 10-k per band)
c  the nir region has three water vapor bands, ten k's for each band.
c    1.   1000-4400 (/cm)         2.27-10.0 (micron)
c    2.   4400-8200               1.22-2.27
c    3.   8200-14300              0.70-1.22
c
c  input parameters:                           units
c    wh,ta,taucl,csm,ibnd,fice,kaer,kprf,idxc,cmix,denn,rh,albb,albd
c    icwp,cwp,cip,cclv,rew,rei,zth
c  fixed input data:
c    h2o absorption coefficient (xk)           cm**2/gm
c    k-distribution function    (hk)           fraction
c
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
cc   99-04-27   yu-tai hou      - updated cloud radiative property
c            calculations use linear t-adjusted method.
c   99-09-13   yu-tai hou      - updated to chou's june,1999 version
c   99-10-13   yu-tai hou      - changed aerosol to opac algorithm
c
c********************************************************************
c
c --- input
                      d i m e n s i o n
     1  wh(imbx,l),   taucl(imbx,l), csm(imax),   rh (imbx,l)
     2, albb(imax),   albd(imax),    zth(imax)
     3, cwp(imbx,l),  cip(imbx,l),   rew(imbx,l), rei(imbx,l)
     4, ccly(imbx,l), ta(imbx,lp1),  fice(imbx,l)
     5, rewi(imbx,l), reii(imbx,l)
     6, idxc(nxc,imax),cmix(nxc,imax),denn(ndn,imax),kprf(imax)
      logical daytm(imax)
c --- output
                      d i m e n s i o n
     1  fnet0 (imbx,lp1), dwsfb0(imax), dwsfd0(imax)
     2, fnetc (imbx,lp1), dwsfbc(imax), dwsfdc(imax)
     3, tupfxc(imax),     supfxc(imax), sdnfxc(imax)
     4, tupfx0(imax),     supfx0(imax), sdnfx0(imax)
!
      integer ncloud
      logical cloudy(imax)
!
c --- temporary array
      real
     1  upflux(imbx,lp1), dwflux(imbx,lp1)
     2, dwsfxb(imax),     dwsfxd(imax)
     3, tauto (imbx,l),   ssato (imbx,l),     asyto (imbx,l)
     4, taurs (l),        ssat1 (imbx,l),     asyt1 (imbx,l)
     5, tauaer(imbx,l),   ssaaer(imbx,l),     asyaer(imbx,l)
     6, xk  (nk0),        hk  (nk0,nrb)
c --- t adjusted cld property method
     6, ssaw0(nrb,2),   ssai0(nrb,2),   asyw0(nrb,2), asyi0(nrb,2)
     7, fffcw(imbx,l),  ffft1(imbx,l),  fffto(imbx,l)
     8, asycw(imbx,l),  ssacw(imbx,l)
                      d i m e n s i o n
     1  a0w(nrb,2), a1w(nrb,2), b0w(nrb,2), b1w(nrb,2)
     2, a0i(nrb,2), a1i(nrb,2), c0w(nrb,2), c1w(nrb,2)
     3, b0i(nrb,2), b1i(nrb,2), b2i(nrb,2), facw(imax,l)
     4, c0i(nrb,2), c1i(nrb,2), c2i(nrb,2), faci(imax,l)
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
c     data ssaw0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9916/
      data ssaw0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9965/
     1,    asyw0/.8678,.8185,.8354,.8315, .8723,.8182,.8354,.8311/
c    2,    ssai0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9750/
     2,    ssai0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9823/
     3,    asyi0/.9058,.8322,.8068,.8220, .9070,.8304,.8067,.8174/
      data fffrs0 / 0.1 /
      data fpmin,fpmax /1.0e-8, 0.999999/
c
                         d a t a
c ---- t-adjusted cld prop coeff, water cloud
     1    a0w / 1.466e-2, 2.276e-2, 2.654e-2, 2.494e-2
     1,         1.528e-2, 2.286e-2, 2.642e-2, 2.517e-2 /
     2,   a1w / 1.617e+0, 1.451e+0, 1.351e+0, 1.392e+0
     2,         1.611e+0, 1.449e+0, 1.353e+0, 1.386e+0 /
     3,   b0w / 1.708e-1, 5.314e-4,-4.594e-6, 6.473e-3
     3,         1.674e-1, 5.427e-4,-3.306e-6, 3.218e-3 /
     4,   b1w / 7.142e-3, 1.258e-3, 2.588e-5, 6.649e-4
     4,         7.561e-3, 1.263e-3, 2.287e-5, 5.217e-4 /
     5,   c0w / 8.266e-1, 7.507e-1, 7.925e-1, 7.811e-1
     5,         8.344e-1, 7.501e-1, 7.922e-1, 7.808e-1 /
     6,   c1w / 4.119e-3, 6.770e-3, 4.297e-3, 5.034e-3
     6,         3.797e-3, 6.812e-3, 4.323e-3, 5.031e-3 /
c
                         d a t a
c --- t-adjusted cld prop coeff, ice cloud
     1    a0i / 2.822e-4,-3.248e-5,-3.758e-5,-1.214e-5
     1,         2.712e-4,-4.308e-5,-3.917e-5,-2.456e-5 /
     2,   a1i / 2.491e+0, 2.522e+0, 2.522e+0, 2.520e00
     2,         2.489e+0, 2.523e+0, 2.522e+0, 2.521e00 /
     3,   b0i / 1.853e-1, 2.544e-3,-7.701e-7, 1.461e-2
     3,         1.738e-1, 2.461e-3,-8.979e-7, 7.083e-3 /
     4,   b1i / 1.841e-3, 1.023e-3, 9.849e-6, 4.612e-4
     4,         1.887e-3, 9.436e-4, 8.102e-6, 3.495e-4 /
     5,   b2i /-6.671e-6,-2.266e-6,-.3988e-9,-1.202e-6
     5,        -6.615e-6,-2.107e-6,-.1862e-9,-8.500e-7 /
     6,   c0i / 8.388e-1, 7.572e-1, 7.519e-1, 7.600e-1
     6,         8.414e-1, 7.566e-1, 7.519e-1, 7.566e-1 /
     7,   c1i / 1.519e-3, 1.563e-3, 1.099e-3, 1.275e-3
     7,         1.477e-3, 1.537e-3, 1.097e-3, 1.241e-3 /
     8,   c2i /-6.702e-6,-5.232e-6,-3.081e-6,-4.020e-6
     8,        -6.403e-6,-5.130e-6,-3.070e-6,-3.804e-6 /
c
      save xk, hk, ssaw0,ssai0, asyw0,asyi0, fffrs0, fpmin, fpmax
      save a0w,a1w,b0w,b1w,c0w,c1w,a0i,a1i,b0i,b1i,b2i,c0i,c1i,c2i
c
      do k=1,l
        do i=1,ipts
          facw(i,k) = max(0.0, min(10.0,273.15-ta(i,k)))*0.1
          faci(i,k) = max(0.0, min(30.0,263.15-ta(i,k)))/30.0
        enddo
      enddo
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
      do ib=ibb1,ibb2
c
c===> ... get aerosols and rayleigh scattering optical properties
c
      ib1 = nvb + ib
      call aeros99(ipts,ib1,kaer,kprf,idxc,cmix,denn,rh
     1,          tauaer,ssaaer,asyaer,taurs)
      cloudy(:) = .false.
c
      if (icwp .eq. 1) then
        do k=1,l
        do i=1,ipts
          if (ccly(i,k) .gt. 0.0) then
c --- t-adj method
            tau1=cwp(i,k)*(  facw(i,k) *(a0w(ib,1)+a1w(ib,1)*rewi(i,k))
     1                 +(1.0-facw(i,k))*(a0w(ib,2)+a1w(ib,2)*rewi(i,k)))
            tau2=cip(i,k)*(  faci(i,k) *(a0i(ib,1)+a1i(ib,1)*reii(i,k))
     1                 +(1.0-faci(i,k))*(a0i(ib,2)+a1i(ib,2)*reii(i,k)))
            taucl(i,k) = tau1 + tau2
            ssa1 = 1.0 - (  facw(i,k) *(b0w(ib,1)+b1w(ib,1)*rew(i,k))
     1               + (1.0-facw(i,k))*(b0w(ib,2)+b1w(ib,2)*rew(i,k)))
            ssa2 = 1.0 - (  faci(i,k) *(b0i(ib,1)
     1                 + (b1i(ib,1)+b2i(ib,1)*rei(i,k))*rei(i,k))
     2               + (1.0-faci(i,k))*(b0i(ib,2)
     3                 + (b1i(ib,2)+b2i(ib,2)*rei(i,k))*rei(i,k)) )
            ssaw1 = ssa1 * tau1
            ssaw2 = ssa2 * tau2
            ssacw(i,k) = ssaw1 + ssaw2
            asy1 =    facw(i,k) *(c0w(ib,1)+c1w(ib,1)*rew(i,k))
     1         + (1.0-facw(i,k))*(c0w(ib,2)+c1w(ib,2)*rew(i,k))
            asy2 =    faci(i,k) *(c0i(ib,1)
     1              + (c1i(ib,1)+c2i(ib,1)*rei(i,k))*rei(i,k))
     2         + (1.0-faci(i,k))*(c0i(ib,2)
     3              + (c1i(ib,2)+c2i(ib,2)*rei(i,k))*rei(i,k))
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
      else
        do k=1,l
        do i=1,ipts
          if (ccly(i,k) .gt. 0.0) then
            tau2 = fice(i,k) * taucl(i,k)
            tau1 = taucl(i,k) - tau2
            ssa1 =      facw(i,k) * ssaw0(ib,1)
     1           + (1.0-facw(i,k))* ssaw0(ib,2)
            ssa2 =      faci(i,k) * ssai0(ib,1)
     1           + (1.0-faci(i,k))* ssai0(ib,2)
            ssaw1 = ssa1 * tau1
            ssaw2 = ssa2 * tau2
            ssacw(i,k) = ssaw1 + ssaw2
            asy1 =      facw(i,k) * asyw0(ib,1)
     1           + (1.0-facw(i,k))* asyw0(ib,2)
            asy2 =      faci(i,k) * asyi0(ib,1)
     1           + (1.0-faci(i,k))* asyi0(ib,2)
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
      end if
!
      ncloud = 0
      do i=1,ipts
        if (cloudy(i)) ncloud = ncloud + 1
      enddo
c
c===> ... ik is the index for the k-distribution function (or the
c     absorption coefficient)
c
        do ik=1,nk0
c
         if (hk(ik,ib) .ge. 0.00001) then
c
c===> ... compute tatal optical thickness, single scattering albedo,
c         and asymmetry factor for clear sky
c
           do k=1,l
           do i=1,ipts
             tauwv      = xk(ik)*wh(i,k)
             tauto(i,k) = max(fpmin, tauwv+tauaer(i,k)+taurs(k))
             ssat1(i,k) = ssaaer(i,k)*tauaer(i,k)+taurs(k)
             asyt1(i,k) = asyaer(i,k)*ssaaer(i,k)*tauaer(i,k)
             ffft1(i,k) = asyaer(i,k)*asyt1(i,k) + fffrs0*taurs(k)
             ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
             tem        = 1.0 / max(fpmin, ssat1(i,k))
             asyto(i,k) = asyt1(i,k) * tem
             fffto(i,k) = ffft1(i,k) * tem
           enddo
           enddo
c
c===> ... clear sky fluxes calculations
c
           call swflux99(ipts,tauto,ssato,asyto,fffto,csm,zth,albb,albd,
     1                 daytm,upflux,dwflux,dwsfxb,dwsfxd)
c
           do k=1,lp1
             do i=1,ipts
               fnet0 (i,k) = fnet0 (i,k)
     1                     + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
             enddo
           enddo
           do i=1,ipts
             tupfx0(i) = tupfx0(i) + upflux(i,1)   * hk(ik,ib)
             supfx0(i) = supfx0(i) + upflux(i,lp1) * hk(ik,ib)
             sdnfx0(i) = sdnfx0(i) + dwflux(i,lp1) * hk(ik,ib)
             dwsfb0(i) = dwsfb0(i) + dwsfxb(i)     * hk(ik,ib)
             dwsfd0(i) = dwsfd0(i) + dwsfxd(i)     * hk(ik,ib)
           enddo
           if (ncloud .gt. 0) then
c
c===> ... compute tatal optical thickness, single scattering albedo,
c         and asymmetry factor for cloudy sky
c
           do k=1,l
             do i=1,ipts
               if (ccly(i,k) .gt. 0.0) then
                 tauto(i,k) = taucl(i,k) + tauto(i,k)
                 ssat1(i,k) = ssacw(i,k) + ssat1(i,k)
                 ssato(i,k) = min(fpmax, ssat1(i,k)/tauto(i,k))
                 tem        = 1.0 / max(fpmin, ssat1(i,k))
                 asyto(i,k) = (asycw(i,k) + asyt1(i,k)) * tem
                 fffto(i,k) = (fffcw(i,k) + ffft1(i,k)) * tem
               end if
             enddo
           enddo
c
c===> ... cloudy sky fluxes calculations
c
           call swflux99(ipts,tauto,ssato,asyto,fffto,csm,zth,albb,albd,
     1                 daytm,upflux,dwflux,dwsfxb,dwsfxd)
c
           do k=1,lp1
             do i=1,ipts
               fnetc(i,k) = fnetc(i,k)
     1                    + (dwflux(i,k) - upflux(i,k))*hk(ik,ib)
             enddo
           enddo
           do i=1,ipts
             tupfxc(i) = tupfxc(i) + upflux(i,1)   * hk(ik,ib)
             supfxc(i) = supfxc(i) + upflux(i,lp1) * hk(ik,ib)
             sdnfxc(i) = sdnfxc(i) + dwflux(i,lp1) * hk(ik,ib)
             dwsfbc(i) = dwsfbc(i) + dwsfxb(i)     * hk(ik,ib)
             dwsfdc(i) = dwsfdc(i) + dwsfxd(i)     * hk(ik,ib)
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
         endif
        enddo           ! k-distribution loop ends here
      enddo           ! loop over nir bands ends here
c
      return
      end
