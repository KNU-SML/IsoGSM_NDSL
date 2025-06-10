#define G g
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
#define PI pi
#include <define.h>
      subroutine rdlwsw(ix,ix2,qs,sigi,sigl,
     1                  lat,xlat,solc,rsin1,rcos1,rcos2,slmsk,
     2                  qqh2o,tt,o3qo3,tsfc,coszro,cldary,cldtot,cldcnv,
#ifndef SW99UPDATE
     3                  albedo,albvb,albnb,albvd,albnd,paers,
#else
     3                  albedo,albvb,albnb,albvd,albnd,
     3                  kprf,idxc,cmix,denn,rhrh,
#endif
     4                  itimsw,itimlw,ko3,kalb,iswsrc,
     4                  ibnd,icwp,qci,icfc,
     5                  htlw,slwup,slwdn,tlwup,
     6                  htsw,sswup,sswdn,tswup,tswdn,
#ifdef CLR
     6                  tlwup0,tswup0,sswdn0,sswup0,slwdn0,
#endif
#ifdef CONDENS
     +                  taucld,cldwp,cldip,
#endif
     7                  sswfvb,sswfnb,sswfvd,sswfnd,solinr)
c
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c
c
c [March 2000] Fanglin Yang
c   add the option to use CHOU_UIUC longwave radiation routine
c -- input preparation
c   taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
c   cld0 (i,l)   - unweighted cloud fraction, k=1-->top layer
c   pl2  (i,lp1) - interfacing pressure (mb), k=1-->model top
c
c-paramodel
c  ********************************************************************
c  *  rewrite the gfdl code "radfs" to call m.d.chous sw radiation
c  *    and to call the gfdl lw radiation         y-t hou dec1994
c  *  updates made to opac aerosol algorithm (1998)...y.h...sep99
c
c  *
c  *  argument list:
c  *    input
c  *      ix,ix2      - first dimensions of arrays
c  *      qs          - surface pressure in cb
c  *      sigi        - model sigma interface level (k=1 at the sfc)
c  *      sigl        - model sigma layer mean value
c  *      lat,xlat    - lat index and radians
c  *      solc        - solar constant in ly/min
c  *      rsin1,rcos1,rcos2
c  *                  - sin and cos lat for ozone interpolations
c  *      slmsk       - land/sea/ice mask (0:sea.1:land,2:ice)
c  *      qqh2o       - sepcific humidity in g/g   (k=1 at the sfc)
c  *      tt          - temperature in k           (k=1 at the sfc)
c  *      o3qo3       - ozone concentration in g/g (k=1 at the sfc)
c  *      tsfc        - surface temperature in k
c  *      coszro      - cosine of zenith angle
c  *      cldary      - packed cloud array         (k=1 at the sfc)
c  *      cldtot      - stratiform cloud
c  *      cldcnv      - convective cloud
c  *      albedo      - surface albedo from climotology
c  *      albvb,albvd - vis band albedoes for beam and diff radiation
c  *      albnb,albnd - nir band albedoes for beam and diff radiation
#ifndef SW99UPDATE
c  *      paers       - aerosol profiles (in fractions)
#else
c  *      kprf        - tropospheric aerosol profile type index
c  *      idxc,cmix   - aerosol component index and mixing ratio
c  *      denn        - aerosol number densities of 1st and 2nd layers
c  *      rhrh        - relative humidity in fraction
#endif
c  *    control flags
c  *      itimsw,itimlw
c  *                  - sw, lw radiation calls duration in hour
c  *      ko3         - ozone data, =1 input data; =0 gfdl climotology
c  *      kalb        - sfc alb, =0 climotology, comp over oceans
c  *                    =1 input four components from calling progm
c  *      iswsrc      - flags for selection of sw absorbers
c  *                    1:aerosols, 2:o2, 3:co2, 4:h2o, 5:o3
c  *                    =0:without; =1 with
c  *    output
c  *      htlw        - lw heating rates in k/sec
c  *      slwup       - sfc upward lw flux in w/m**2
c  *      slwdn       - sfc downward lw flux in w/m**2
c  *      tlwup       - toa upward lw flux in w/m**2
c  *      htsw        - sw heating rates in k/sec
c  *      sswup       - sfc upward sw flux in w/m**2
c  *      sswdn       - sfc downward sw flux in w/m**2
c  *      tswup       - toa upward sw flux in w/m**2
c  *      tswdn       - toa downward sw in w/m**2
c  *      tlwup0      - clear sky toa upward lw flux
c         tswup0      - clear sky toa upward sw flux
c  *      sswdn0      - clear sky sfc downward sw flux
c  *      sswup0      - clear sky sfc upward sw flux
c  *      slwdn0      - clear sky sfc downward lw flux
c  *      sswfvb      - vis beam down sw flux at sfc in w/m**2
c  *      sswfvd      - vis diff down sw flux at sfc in w/m**2
c  *      sswfnb      - nir beam down sw flux at sfc in w/m**2
c  *      sswfnd      - nir diff down sw flux at sfc in w/m**2
c  ********************************************************************
c
c prgmmr:  yifeng cui add rsm mp     org:sdsc     date:02-09-09
c
#ifndef SW99UPDATE
#ifndef ICECLOUD
#include <rdparm.h>
#else
#include <rdparm8.h>
#endif
#else
#include <rdparm99.h>
#include <aerparm.h>
#endif
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define LONF2S igrd12p_*jgrd12p_
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONF2S igrd12p_
#define LONLENS lonlen(mype)*2
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
      parameter (lngth=37*l)
c --- input from from calling program
                      d i m e n s i o n
     1  qs   (ix2),  sigi (ix2,lp1),  sigl (ix2,l),    slmsk(ix2)
     2, qqh2o(ix2,l),tt   (ix2,l),tsfc (ix2),   o3qo3(ix2,l)
     3, cldary(ix2,l),cldtot(ix2,l),cldcnv(ix2,l)
     3, coszro(ix2),albedo(ix2),xlat(ix2)
     4, albvb(ix2),   albnb(ix2),   albvd(ix2),   albnd(ix2)
#ifndef SW99UPDATE
     5, paers(ix2,nae-1),          iswsrc(nsrc),  qci(ix2,l)
#else
     5, rhrh(ix2,l), iswsrc(nsrc), qci(ix2,l)
     6, denn(ndn,ix2),idxc(nxc,ix2),cmix(nxc,ix2),kprf(ix2)
#endif
c --- output to calling program
                      d i m e n s i o n
     1  htlw (ix2,l), slwup(ix2), slwdn(ix2), tlwup(ix2)
     2, htsw (ix2,l), sswup(ix2), sswdn(ix2), tswup(ix2), tswdn(ix2)
c --- four components of downward sw flux
     3, sswfvb(ix2), sswfnb(ix2),sswfvd(ix2),sswfnd(ix2)
#ifdef CLR
     4, tlwup0(ix2), tswup0(ix2),sswdn0(ix2),sswup0(ix2),slwdn0(ix2)
#endif
     5, solinr(ix2,l)
#ifdef CONDENS
      dimension  taucld(ix2,l),cldwp(ix2,l),cldip(ix2,l)
#endif
c
c --- internal arrays
                      d i m e n s i o n
c --- local arrays for radiative quantities
     1  hlw(imbx,l), tlwuc(imax), slwnc(imax), tswdc(imax)
     2, hsw(imbx,l), tswuc(imax), sswuc(imax), sswdc(imax)
     3, tlwu0(imax), tswu0(imax), sswd0(imax), sswu0(imax), slwn0(imax)
c --- local arrays for clouds
                      d i m e n s i o n
     1  nclds(imax),     cfacsw(imbx,lp1), cfaclw(imbx,lp1,lp1)
     2, ktop (imbx,lp1), kbtm  (imbx,lp1)
     3, cldlw(imbx,lp1), cldsw (imbx,l),   taucl (imbx,l)
     4, cwp  (imbx,l),   cip   (imbx,l),   fice  (imbx,l)
     5, rew  (imbx,l),   rei   (imbx,l)
     6, taulw(imbx,l),   cld0(imbx,l) , clwp(imbx,l), solinsw(imbx,l)
c --- local arrays for others
                      d i m e n s i o n
     1  prssi(imbx,lp1), prssl(imbx,lp1), temp (imbx,lp1),sfcalb(imax)
#ifndef SW99UPDATE
     2, paera(imbx,nae), rh2o (imbx,l),   qo3  (imbx,l),  coszen(imax)
#else
     2, rhrh1(imbx,l),   rh2o (imbx,l),   qo3  (imbx,l),  coszen(imax)
     2, denn1(ndn,imax), idxc1(nxc,imax), cmix1(nxc,imax),kprf1(imax)
#endif
     3, alvb (imax),     alnb (imax),     alvd (imax),    alnd (imax)
     4, gdfvb(imax),     gdfnb(imax),     gdfvd(imax),    gdfnd(imax)
     5, jjrow(imax),     do3v (imax),     do3vp(imax),    tthan(imax)
     6, pl2  (imbx,lp1)
c --- seasonal climatologies of o3 (obtained from a previously run
c        code which interpolates o3 to user vertical coordinate).
c        defined as 5 deg lat means n.p.->s.p.
                      c o m m o n /savmem/
c --- ...winter....  ...spring....  ...summer....  ....fall.....
     1  dduo3n(37,l), ddo3n2(37,l), ddo3n3(37,l), ddo3n4(37,l)
                      d i m e n s i o n
     1  rad1(lngth),  rad2(lngth),  rad3(lngth),  rad4(lngth)
      equivalence (rad1(1),dduo3n(1,1)),(rad2(1),ddo3n2(1,1))
      equivalence (rad3(1),ddo3n3(1,1)),(rad4(1),ddo3n4(1,1))
c --- surface albedo (used for the old scheme)
                      c o m m o n /ssalb/
     1  albd(21,20),  za(20),       trn(21),      dza(19)
      common /rdfsav/  degrad,hsigma,daysec,rco2
c===> ... begin here
c       solc,the solar constant is scaled to a more current value.
c       i.e. if solc=2.0 ly/min then ssolar=1.96 ly/min. then
c       convert unit to w/m**2
      ssolar = solc * 0.98e0 * 6.97667e2
c
      nrchnk=(ix-1)/imax+1
      do 4000 irchnk=1,nrchnk
        ibeg = (irchnk - 1) * imax + 1
        iend = ibeg + imax - 1
        if (iend .gt. ix) iend = ix
        ipts = iend - ibeg + 1
c
c===> ... assign temp,pressures, cosin of zenith angle
c         note: the nmc variables are in mks, gfdl lw variables
c         are in cgs units, so prssl is in dyns/cm**2. but
c         sw prssi is in mb.
      do 20 i=1,ipts
        ir = i + ibeg - 1
        prssl(i,lp1) = 1.0e4 * qs(ir)
        prssi(i,lp1) = 10.0e0 * qs(ir)
        pl2  (i,lp1) = 10.0e0 * qs(ir)
        temp (i,lp1) = max(100.,tsfc(ir))
        coszen(i)    = coszro(ir)
#ifndef SW99UPDATE
c===> ... currently the convective type of aerosol is not in use
c         to give zeroes here.
        paera(i,6)   = 0.0e0
#else
        kprf1(i)     = kprf(ir)
#endif
  20  continue
#ifndef SW99UPDATE
      do 30 k=1,5
      do 30 i=1,ipts
        ir = i + ibeg - 1
        paera(i,k) = paers(ir,k)
  30  continue
#else
      do 30 i=1,ipts
        ir = i + ibeg - 1
        do j=1,nxc
          idxc1(j,i) = idxc(j,ir)
          cmix1(j,i) = cmix(j,ir)
        end do
        do j=1,ndn
          denn1(j,i) = denn(j,ir)
        end do
  30  continue
#endif
c===> ... all radiation variables have k=1 at the top of the atmosphere.
c         while the forcast model variables have k=1 at the surface.
      do 50 k=1,l
      do 50 i=1,ipts
        ir = i + ibeg - 1
        temp (i,k) = max(100.,tt(ir,lp1-k))
        prssl(i,k) = 1.0e4 * sigl(i,lp1-k) * qs(ir)
        prssi(i,k) = 10.0e0 * sigi(i,lp2-k) * qs(ir)
        pl2  (i,k) = 10.0e0 * sigi(i,lp2-k) * qs(ir)
        rh2o (i,k) = max(3.0e-6, qqh2o(ir,lp1-k))
#ifdef SW99UPDATE
        rhrh1(i,k) = rhrh(ir,lp1-k)
#endif
  50  continue
      if(icwp.eq.1) then
       tmp1 = 1.0e5 / g_
        do k=1,l
          do i=1,ipts
            ir = i + ibeg - 1
            rhodz = (prssi(i,k+1)-prssi(i,k))*tmp1
            clwp (i,k) = max(qci(ir,lp1-k),0.) * rhodz          ! g/m**2
c-------------------------
c for t-inversion-type cloud 
c clwp for t-inversion-type cloud is from Duynkerke and Teixeira
c (2001,JC,p1466-1478),monthly mean value in table1 
c-------------------------
            icl=nint(cldary(ir,lp1-k))
            if (icl.eq.2) then
               clwp (i,k) = max(clwp(i,k),57.0e0)         ! g/m**2
            endif
c-------------------------
c for convective cloud
c clwp for (deep) convective cloud is from...
c-------------------------
            if (cldcnv(ir,lp1-k).gt.0.0) then
              clwp(i,k)=30.0   ! g/m**2
            endif
          enddo
        enddo
      endif
c****************************************************
c         ozone section
c****************************************************
      if (ko3 .gt. 0) then
        do 60 k=1,l
        do 60 i=1,ipts
c ken+moorthi: set up minimum positive value for o3
          qo3(i,k) = max(o3qo3(i+ibeg-1,lp1-k), 1.0e-10)
c         qo3(i,k) = o3qo3(i+ibeg-1,lp1-k)
  60    continue
      else
        do 80 i=1,ipts
          ir = i + ibeg - 1
          th2= 0.2e0 * xlat(ir) * degrad
          jjrow(i) = 19.001e0 - th2
          tthan(i) = (19-jjrow(i)) - th2
  80    continue
c===> ... compute climatological zonal mean ozone,
c         seasonal and spatial interpolation done below.
        do 100 k=1,l
        do 100 i=1,ipts
          do3v(i)  = dduo3n(jjrow(i),k) + rsin1*ddo3n2(jjrow(i),k)
     1             + rcos1*ddo3n3(jjrow(i),k)
     2             + rcos2*ddo3n4(jjrow(i),k)
          do3vp(i) = dduo3n(jjrow(i)+1,k) + rsin1*ddo3n2(jjrow(i)+1,k)
     1             + rcos1*ddo3n3(jjrow(i)+1,k)
     2             + rcos2*ddo3n4(jjrow(i)+1,k)
c===> ... now latitudinal interpolation, and
c         convert o3 into mass mixing ratio(original data mpy by 1.e4)
          qo3(i,k) = 1.0e-4 * (do3v(i)+tthan(i)*(do3vp(i)-do3v(i)))
 100    continue
      end if
c******************************************************
c         cloud optical properties section
c******************************************************
c
#ifdef CONDENS
       do k=1,l
       do i=1,ipts
         ir = i + ibeg - 1
         taucld(ir,k) = 0.0
         cldwp(ir,k)  = 0.0
         cldip(ir,k)  = 0.0
       enddo
       enddo
#endif
c
cc     call cldprp(ix2,prssi,temp,cldary,ibeg,ipts,xlat,
cc   1             ktop,kbtm,nclds,cldlw,taucl,cfacsw,cldsw)
#ifndef SW99UPDATE
#ifndef ICECLOUD
       call cldprp(ix2,prssi,temp,cldary,ibeg,ipts,xlat,
     1             ktop,kbtm,nclds,cldlw,taucl,cfacsw,cldsw
     2            ,taulw,cld0)
#else
c      change the unit of prssl to [hPa] same as prssi
       prssl(1:ipts,1:lp1)=prssl(1:ipts,1:lp1)*1.e-3
       call cldprpi
     *         (ix2,sigl,prssl,prssi,temp,cldtot,cldcnv,ibeg,ipts,lat
     1,         icwp,clwp,slmsk,xlat,ktop,kbtm,nclds
     2,         cldlw,taucl,cfacsw,cldsw
     3,         cwp,cip,rew,rei,fice
     4,         taulw,cld0)
#endif
#else
       call cldprp00(ix2,prssi,temp,cldtot,cldcnv,ibeg,ipts,lat
     1,            icwp,clwp,slmsk,xlat,ktop,kbtm,nclds
     2,            cldlw,taucl,cfacsw,cldsw
     3,            cwp,cip,rew,rei,fice
     4,            taulw,cld0)
#endif
c
#ifdef CONDENS
       do k=1,l
       do i=1,ipts
         ir = i + ibeg - 1
         kk = l - k + 1                 ! re-ordering (sfc->top)
         taucld(ir,kk) = taulw(i,k)     ! tau_0
         if (icwp.eq.1) then
          cldwp(ir,kk) = cwp(i,k)*0.001 ! kg/m**2
          cldip(ir,kk) = cip(i,k)*0.001 ! kg/m**2
         endif
       enddo
       enddo
#endif
c
c******************************************************
c         surface albedo section
c******************************************************
c
      if (kalb .eq. 0) then
        iq = int(20.0e0 * 0.537e0 + 1.0e0)
        do 110 i=1,ipts
          ir = i + ibeg -1
          sfcalb(i) = albedo(ir)
c===> ... the following code gets albedo from payne,1972 tables if
c      1) open sea point (slmsk=0);   2) kalb=0
          if(coszen(i).gt.0.0 .and. slmsk(ir).eq.0.0) then
            zen = degrad*acos(max(coszen(i),0.0e0))
            if(zen .ge. 74.e0) jx = int(0.5e0*(90.e0-zen)+1.e0)
            if(zen .lt. 74.e0 .and. zen .ge. 50.e0)
     1        jx = int(0.25e0*(74.e0-zen) + 9.0e0)
            if(zen.lt.50.e0) jx = int(.1e0*(50.e0-zen)+15.0e0)
            dzen = -(zen-za(jx))/dza(jx)
            alb1 = albd(iq,jx)+dzen*(albd(iq,jx+1)-albd(iq,jx))
            alb2 = albd(iq+1,jx)+dzen*(albd(iq+1,jx+1)-albd(iq+1,jx))
            sfcalb(i) = alb1+20.e0*(alb2-alb1)*(0.537e0-trn(iq))
          endif
 110    continue
        do 120 i=1,ipts
          alvd(i) = sfcalb(i)
          alnd(i) = sfcalb(i)
          alvb(i) = sfcalb(i)
          alnb(i) = sfcalb(i)
 120    continue
c===> ... visible and near ir direct beam albedo,if not ocean nor snow
c         function of cosine solar zenith angle..
        do 130 i=1,ipts
          if (slmsk(i+ibeg-1).gt.0.0e0 .and. sfcalb(i).le.0.5e0) then
            albd0 = -18.0e0 * (0.5e0 - acos(coszen(i))/pi_)
            albd0 = exp (albd0)
            alvd1 = (alvd(i) - 0.054313e0) / 0.945687e0
            alvb(i) = alvd1 + (1.0e0 - alvd1) * albd0
            alnb(i) = alvb(i)
          end if
 130    continue
      else
        do 150 i=1,ipts
          ir = i + ibeg - 1
          alvd(i) = albvd(ir)
          alnd(i) = albnd(ir)
          alvb(i) = albvb(ir)
          alnb(i) = albnb(ir)
 150    continue
      end if
c***********************************************
c         final check before radiation calls
c***********************************************
      if (ipts .lt. imax) then
        ipts1 = ipts + 1
        do 160 i=ipts1,imax
          coszen(i) = coszen(ipts)
          alvd(i) = alvd(ipts)
          alnd(i) = alnd(ipts)
          alvb(i) = alvb(ipts)
          alnb(i) = alnb(ipts)
 160    continue
        do 170 k=1,lp1
        do 170 i=ipts1,imax
          prssl(i,k) = prssl(ipts,k)
          prssi(i,k) = prssi(ipts,k)
          pl2  (i,k) = pl2  (ipts,k)
          temp (i,k) = temp (ipts,k)
 170    continue
        do 180 k=1,l
        do 180 i=ipts1,imax
          rh2o(i,k) = rh2o(ipts,k)
          qo3 (i,k) = qo3 (ipts,k)
 180    continue
      end if
c
      if(itimsw .eq. 0) go to 300
c*************************************************
c         calling chous sw radiation routine
c*************************************************
c
c===> ... set timer for sw rad call
c
c     t00  = second()
c     ovhd = second() - t00
c     tbef = second()
c
#ifndef SW99UPDATE
      call swr95(ipts,ssolar,iswsrc,ibnd,prssi,
     1           temp,rh2o,qo3,rco2,coszen,taucl,
     2           cldsw,cfacsw,icfc,icwp,
#ifdef ICECLOUD
     2           cwp,cip,rew,rei,fice,
#endif
     3           alvb,alvd,alnb,alnd,paera,
     4           hsw,tswuc,tswdc,sswuc,sswdc,
     5           tswu0,sswu0,sswd0,
     6           gdfvb,gdfvd,gdfnb,gdfnd,solinsw)
#else
      call swr99(ipts,ssolar,iswsrc,ibnd,prssi,
     1           temp,rh2o,qo3,rco2,coszen,taucl,
     2           cldsw,cfacsw,icfc,icwp,cwp,cip,rew,rei,fice,
     3           alvb,alvd,alnb,alnd,kprf1,idxc1,cmix1,denn1,rhrh1,
     4           hsw,tswuc,tswdc,sswuc,sswdc,
     5           tswu0,sswu0,sswd0,
     6           gdfvb,gdfvd,gdfnb,gdfnd)
#endif
c
c===> ... save toa and sfc fluxes in w/m**2
      do 220 i=1,ipts
        ir = i + ibeg - 1
        tswup (ir) = tswuc(i)
        tswdn (ir) = tswdc(i)
        sswup (ir) = sswuc(i)
        sswdn (ir) = sswdc(i)
#ifdef CLR
        tswup0(ir) = tswu0(i)
        sswdn0(ir) = sswd0(i)
        sswup0(ir) = sswu0(i)
#endif
c===> ... diffused downward sfc fluxes (vis,nir)
        sswfvd(ir) = gdfvd(i)
        sswfnd(ir) = gdfnd(i)
c===> ... direct beam downward sfc fluxes (vis,nir)
        sswfvb(ir) = gdfvb(i)
        sswfnb(ir) = gdfnb(i)
 220  continue
c===> ... convert heating rates to deg/sec
      do 240 k=1,l
      do 240 i=1,ipts
        htsw(i+ibeg-1,lp1-k) = hsw(i,k) * daysec
 240  continue
c===> ... tau (vis), solinr
      do k=1,l
      do i=1,ipts
cc        tau(i+ibeg-1,lp1-k) = taucl(i,k)
        solinr(i+ibeg-1,lp1-k) = solinsw(i,k)
      end do
      end do
                                                                                
c
 300  if(itimlw .eq. 0) go to 400
#ifndef RADMDC
cc*************************************************
cc        calling gfdl long wave radiation routine
cc*************************************************
cc==> ...  get cld factor for lw calculations
cc
cc    tbef = second()
cc
      call clo89(ipts,cfaclw,cldlw,nclds,kbtm,ktop)
cc
      call lwr88(ipts,hlw,slwnc,tlwuc,
#ifdef CLR
     1               slwn0,tlwu0,
#endif
     1           prssl,temp,rh2o,qo3,cfaclw,cldlw,nclds,ktop,kbtm)
#else
c
c-------------------------------------------------------------
c-------------------------------------------------------------
c     calling uiuc_chou long wave radiation routine

      call lwrchou(ipts,hlw,slwnc,tlwuc,
#ifdef CLR
     1              slwn0,tlwu0,
#endif
     1              pl2,temp,rh2o,qo3,taulw,cld0)
#endif
c-------------------------------------------------------------
c-------------------------------------------------------------
c
c===> ... save toa and sfc fluxes in w/m**2
      do 320 i=1,ipts
        ir = i + ibeg - 1
        tlwup (ir) = 1.0e-3*tlwuc(i)
        slwup (ir) = 1.0e-3*hsigma*temp(i,lp1)**4
        slwdn (ir) = slwup(ir) - 1.0e-3*slwnc(i)
#ifdef CLR
        tlwup0(ir) = 1.0e-3*tlwu0(i)
        slwdn0(ir) = slwup(ir) - 1.0e-3*slwn0(i)
#endif
 320  continue
c===> ... convert heating rates to deg/sec
      do 340 k=1,l
      do 340 i=1,ipts
        htlw(i+ibeg-1,lp1-k) = hlw(i,k) * daysec
 340  continue
 400  continue
c
4000  continue
c
      return
      end
