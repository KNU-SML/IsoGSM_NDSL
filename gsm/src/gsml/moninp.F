#define CP cp
#define G g
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
#define RD rd
#define RV rv
#include <define.h>
      subroutine moninp(ims2,imx2,kmx,ndiff,a,b,tau,rtg,
     1     u1,v1,t1,q1,
     2     pstar,rbsoil,cd,ch,fm,fh,tsea,qss,dphi,spd1,
     3     si,del,sl,slk,rcl,deltim,lat,kdt,thour,
     4     dusfc,dvsfc,dtsfc,dqsfc,hpbl,hgamt,hgamq,slmask,
#ifdef ISOTOPE
     5     revap,dqsfci1,dqsfci2,qssi,tskin,
     $     rresv1,rresv2,xlon,xlat,
#endif
#ifdef CO2PROG
     *     co2fl,dqsfci1,
#endif
#ifndef VIC
     5     kpbl)
#else
     5     kpbl,heat,evap)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    moninp      computes vertical diffusion terms.
c   prgmmrs: mrb personnel   org: w/nmc23    date: 88-04-29
c
c abstract: computes the parameterized effectsc of vertical turbulent
c   eddy diffusion of momentum water vapor and sensible heat in the
c   mrf model. the scheme considers the nonlocal-k approach for the
c   mixed layer and the local-k aproach for the free atmosphere.
c   the mixed layer diffusion is formulated based on the large eddy simulation
c   research. in this approch, the turbulent diffusivity coefficients
c   are calculated from a prescribed profile shape as a function of
c   boundary layer heights and scale parameters considering the
c   countergradient turbulent mixing. in the free atmospheric diffusion,
c   local-k approach is formulated based on scale parameters obtained
c   from aircraft observations.
c
c reference:
c   hong and pan, 1996: nonlocal boundary layer vertical diffusion
c       in a medium-range forecast model. mon. wea. rev., 124, 2322-2339.
c
c program history log:
c   88-04-29  hua-lu pan
c   88-10-28  sela
c   92-09-01  iredell
c   95-09-01  song-you hong for non-local pbl
c   00-04-01  song-you hong for multiple matrix solver and tracer diffusion
c   02-09-09  yifeng cui for rsm mp
c
c usage:    call moninp(im,im2,km,a,b,tau,rtg,
c    1           u1,v1,t1,q1,
c    2           pstar,rbsoil,cd,ch,fm,fh,tsea,qss,dphi,spd1,
c    3           si,del,sl,slk,rcl,deltim,tsea,lat,kdt,thour,
c    4           dusfc,dvsfc,dtsfc,dqsfc,hpbl,hgamt,hgamq)
c   input argument list:
c     im       - number of profiles to compute
c     im2      - first dimension of field slices
c     km       - number of vertical levels
c     u1       - (im2,km) zonal wind * cos(lat) in m/s
c     v1       - (im2,km) merid wind * cos(lat) in m/s
c     t1       - (im2,km) temperature in k
c     q1       - (im2,km) specific humidity in kg/kg
c     ndiff    - number of passive variables including moisture
c     pstar    - (im) surface pressure in kpa
c     rbsoil   - (im) bulk richardson number at surface layer
c     cd       - (im) k**2/(fm*fm) momentum exchange coefficient
c     ch       - (im) k**2/(fm*fh) heat & moisture exchange coefficient
c     fm       - (im) integral of profile function for momentum
c     fh       - (im) integral of profile function for heat & moisture
c     tsea     - (im) surface temperature in k
c     qss      - (im) specific humidity at surface in kg/kg
c     dphi     - (im) coefficient modulating sfc evaporation
c     spd1     - (im) wind speed at lowest sigma level in m/s
c     si       - (km+1) p/psfc at base of layer
c     del      - (km) positive increment of p/psfc across layer
c     sl       - (km) p/psfc at middle of layer
c     slk      - (km) sigma values to the kappa
c     rcl      - reciprocal of square of cos(lat)
c     deltim   - time step in secs
c     lat      - latitude number
c     kdt      - timestep number
c     thour    - forecast hour
c   output argument list:
c     a        - (im,km) negative tendency for v wind in m/s/s
c     b        - (im,km) tendency for u wind in m/s/s
c     tau      - (im,km) tendency for temperature in k/s
c     rtg      - (im,km) tendency for specific humidity in kg/kg/s
c     dusfc    - (im) zonal stress on surface in n/m**2
c     dvsfc    - (im) merid stress on surface in n/m**2
c     dtsfc    - (im) sensible heat flux on surface in w/m**2
c     dqsfc    - (im) latent heat flux on surface in w/m**2
c     hpbl     - (im) pbl height in m
c     hgamt    - (im) countergradient term for temperature in k
c     hgamq    - (im) countergradient term for moisture in kg/kg
c
c subprograms called:
c     tridi2   - solve tridiagonal matrix equation
c
c subprograms called:
c     tridin   - solve tridiagonal matrix equation for t, q, and tracers
c     tridi2   - solve tridiagonal matrix equation for momentum
c
c$$$
#include <paramodel.h>
      parameter(cp=cp_,g=g_,rd=rd_,rv=rv_,hvap=hvap_)
      parameter(gor=g/rd,gocp=g/cp,rocp=rd/cp,fv=rv/rd-1)
      parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g)
      parameter(rlam=150.,vk=0.4,vk2=vk*vk,prmin=1.0,prmax=4.)
      parameter(dw2min=0.0001,dkmin=1.0,dkmax=1000.,rimin=-100.)
      parameter(rbcr=0.5,cfac=7.8,pfac=2.0,sfcfrac=0.1)
      parameter(qmin=1.e-8,xkzo=1.,zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(gamcrt=3.,gamcrq=2.e-3)
      parameter(iun=84)
      dimension dusfc(imx2),dvsfc(imx2),dtsfc(imx2),dqsfc(imx2)
      dimension a(imx2,kmx),b(imx2,kmx),tau(imx2,kmx),
     1     rtg(imx2,kmx*ndiff),
     1     u1(imx2,kmx),v1(imx2,kmx),t1(imx2,kmx),q1(imx2,kmx*ndiff),
     2     tsea(imx2),qss(imx2),pstar(imx2),dphi(imx2),
     3     si(imx2,kmx+1),del(imx2,kmx),sl(imx2,kmx),
     $     slk(imx2,kmx),hpbl(imx2),
     4     fm(imx2),fh(imx2),rbsoil(imx2),cd(imx2),ch(imx2),
     5     hgamt(imx2),hgamq(imx2),slmask(imx2),kpbl(imx2)
#ifdef CO2PROG
** co2 **
      dimension
     $     co2fl(imx2),dqsfci1(imx2)
#endif
c local array
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
      dimension dzot(ILOTS,levs_-1),rdzt(ILOTS,levs_-1),
     $          rdzt1(ILOTS),
     1          betaw(ILOTS),betaq(ILOTS),
     *          betat(ILOTS),
     1          zi(ILOTS,levs_+1),
     *          zl(ILOTS,levs_),
     *          zl1(ILOTS),
     2          dku(ILOTS,levs_-1),
     *          dkt(ILOTS,levs_-1),
     3          al(ILOTS,levs_-1),
     *          ad(ILOTS,levs_),
     *          au(ILOTS,levs_-1),
     4          a1(ILOTS,levs_),
     *          a2(ILOTS,levs_*ndiff),
     5          wscale(ILOTS),
     6          ustar(ILOTS),
     *          spd1(ILOTS),
     *          the1v(ILOTS),
     6          thermal(ILOTS),
     *          rbdn(ILOTS),
     *          rbup(ILOTS),
     7          heat(ILOTS),
     *          evap(ILOTS),
     *          thesv(ILOTS),
     *          the1(ILOTS),
     7          wstar(ILOTS),
#ifndef VIC
     8          phim(ILOTS),
     *          phih(ILOTS)
#else
     8          phim(ILOTS),
     *          phih(ILOTS),
     *          beta2(ILOTS)
#endif
      logical   pblflg(ILOTS),
     *          sfcflg(ILOTS),
     *          stable(ILOTS)
c
#ifdef ISOTOPE
      dimension dqsfci1(imx2),dqsfci2(imx2)
      dimension hgamqi1(imx2),hgamqi2(imx2)
      dimension evapi1(imx2),evapi2(imx2)
      dimension rresv1(ILOTS),rresv2(ILOTS)
      dimension qssi(ILOTS,2)
      dimension revap(ILOTS,2)
      dimension tskin(ILOTS)
      dimension rsea(imx2,2)
      dimension xlon(imx2)
      dimension xlat(imx2)
      real      freq
      real      frkin
      integer   flg(imx2),flg2(imx2)
#endif
c
      im=ims2
      km=kmx
c-----------------------------------------------------------------------
c
 601  format(1x,' moninp lat lon step hour ',3i6,f6.1)
 602      format(1x,'    k','        z','        t','       th',
     1     '      tvh','        q','        u','        v',
     2     '       sp')
 603      format(1x,i5,8f9.1)
 604      format(1x,'  sfc',9x,f9.1,18x,f9.1)
 605      format(1x,'    k      zl    spd2   thekv   the1v'
     1         ,' thermal    rbup')
 606      format(1x,i5,6f8.2)
 607      format(1x,' kpbl    hpbl      fm      fh   hgamt',
     1         '   hgamq      ws   ustar      cd      ch')
 608      format(1x,i5,9f8.2)
 609      format(1x,' k pr dkt dku ',i5,3f8.2)
 610      format(1x,' k pr dkt dku ',i5,3f8.2,' l2 ri t2',
     1         ' sr2  ',2f8.2,2e10.2)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     compute preliminary variables
c
ccc   iprt = 0
ccc   if(iprt.eq.1) then
ccc   latd = 0
ccc   lond = 0
ccc   else
ccc   latd = 0
ccc   lond = 0
ccc   endif
c
      dt    = 4. * deltim
      rdt   = 1. / dt
      do i=1,im
      rdzt1(i) = gor * sl(i,1) / del(i,1)
      enddo
      kmpbl = km / 2
c
      do i=1,im
      do k = 1,km-1
         rdzt(i,k) = gor * si(i,k+1) / (sl(i,k) - sl(i,k+1))
         dzot(i,k) = log(si(i,k+1) / si(i,k)) / gor
      enddo
      enddo
c
      do i = 1,im
         zi(i,1) = 0.
      enddo
      do k = 1, km - 1
         do i = 1, im
            zi(i,k+1) = zi(i,k) - t1(i,k) * dzot(i,k)
         enddo
      enddo
c
      do i = 1,im
         dusfc(i) = 0.
         dvsfc(i) = 0.
         dtsfc(i) = 0.
         dqsfc(i) = 0.
#ifdef ISOTOPE
         dqsfci1(i) = 0.
         dqsfci2(i) = 0.
#endif
#ifdef CO2PROG
         dqsfci1(i) = 0.
#endif
         hgamt(i) = 0.
         hgamq(i) = 0.
         wscale(i) = 0.
         kpbl(i) = 1
         hpbl(i) = zi(i,2)
         pblflg(i) = .true.
         sfcflg(i) = .true.
         if(rbsoil(i).gt.0.0) sfcflg(i) = .false.
      enddo
c
      do i = 1,im
         bet1 = dt*rdzt1(i)*spd1(i)/t1(i,1)
#ifdef VIC
         beta2(i) = dt*rdzt1(i)/t1(i,1)
#endif
         betaw(i) = bet1*cd(i)
         betat(i) = bet1*ch(i)
         betaq(i) = dphi(i)*betat(i)
      enddo
c
      do i = 1,im
         zl1(i) = 0.-(t1(i,1)+tsea(i))/2.*log(sl(i,1))/gor
         ustar(i) = sqrt(cd(i)*spd1(i)**2)
      enddo
c
      do i = 1,im
         thesv(i)  = tsea(i)*(1.+fv*max(qss(i),qmin))
         the1(i)   = t1(i,1)/slk(i,1)
         the1v(i) = the1(i)*(1.+fv*max(q1(i,1),qmin))
         thermal(i) = the1v(i)
         dthe1  = (the1(i)-tsea(i))
#ifdef HYDRO
         dq1    = (q1(i,1) - qss(i))
#else
         dq1    = (max(q1(i,1),qmin) - max(qss(i),qmin))
#endif
#ifndef VIC
         heat(i) = -ch(i)*spd1(i)*dthe1
         evap(i) = -ch(i)*spd1(i)*dq1
#ifdef ISOTOPE
#ifdef NOFRAC
         rsea(i,1)=1.D0
         rsea(i,2)=0.D0
#else
         rsea(i,1)=1.D0
         rsea(i,2)=1.D0
#endif
         if (dq1.le.0.D0) then !! evaporation
#ifdef NOFRAC
         if (xlon(i)/3.1416*180.gt.180 .and.
     $       xlon(i)/3.1416*180.lt.240 .and.
     $       xlat(i)/3.1416*180.gt.10  .and.
     $       xlat(i)/3.1416*180.lt.30 ) then
           revap(i,1)=1.D0 !! sea
           revap(i,2)=0.D0 !! land
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
           flg(i)=2
         elseif (xlon(i)/3.1416*180.gt.220 .and.
     $       xlon(i)/3.1416*180.lt.250 .and.
     $       xlat(i)/3.1416*180.gt.30  .and.
     $       xlat(i)/3.1416*180.lt.50 ) then
           revap(i,1)=0.D0 !! sea
           revap(i,2)=1.D0 !! land
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
           flg(i)=3
         else
           revap(i,1)=0.D0 !! sea
           revap(i,2)=0.D0 !! land
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
         endif
#else
         if (slmask(i).eq.0.D0) then
           qssi(i,1)=(1./freq(1,tskin(i)))*rsea(i,1)*qss(i)
           qssi(i,2)=(1./freq(2,tskin(i)))*rsea(i,2)*qss(i)
           if ((q1(i,1+2*km).gt.0.D0).and.
     $        (q1(i,1).gt.0.D0).and.
     $        (qssi(i,1).gt.0.D0).and.
     $        (qss(i).gt.0.D0).and.
     $        (qssi(i,2).gt.0.D0).and.
     $        (q1(i,1+3*km).gt.0.D0)) then            
           revap(i,1)=frkin(1,spd1(i))
           revap(i,2)=frkin(2,spd1(i))
           dq1    = (max(q1(i,1+2*km),qmin*q1(i,1+2*km)/q1(i,1))
     $          - max(qssi(i,1),qmin*qssi(i,1)/qss(i)))
           evapi1(i)= -ch(i)*spd1(i)*dq1*revap(i,1)
           dq1    = (max(q1(i,1+3*km),qmin*q1(i,1+3*km)/q1(i,1))
     $          - max(qssi(i,2),qmin*qssi(i,2)/qss(i)))
           evapi2(i)= -ch(i)*spd1(i)*dq1*revap(i,2)
           flg(i)=1
           else
           revap(i,1)=1.D0/freq(1,tskin(i))*rsea(i,1)
           revap(i,2)=1.D0/freq(2,tskin(i))*rsea(i,2)
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
           flg(i)=2
           endif
         elseif (slmask(i).ge.1.D0) then !! land and ice
           revap(i,1)=rresv1(i)
           revap(i,2)=rresv2(i)
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
           flg(i)=3
         endif
#endif
         else !! dew
           if (q1(i,1).gt.0.D0) then
           revap(i,1)=freq(1,tskin(i))*q1(i,1+2*km)/q1(i,1)
           revap(i,2)=freq(2,tskin(i))*q1(i,1+3*km)/q1(i,1)
           evapi1(i)=revap(i,1)*evap(i)
           evapi2(i)=revap(i,2)*evap(i)
           flg(i)=4
           else
           revap(i,1)=0.D0
           revap(i,2)=0.D0
           evapi1(i)=0.D0
           evapi2(i)=0.D0
           flg(i)=5
           endif
         endif
#endif
#endif
      enddo
c
c
c     compute the first guess of pbl height
c
      do i = 1, im
         stable(i) = .false.
         zl(i,1) = zl1(i)
         rbup(i) = rbsoil(i)
      enddo
      do k = 2, kmpbl
         do i = 1, im
            if(.not.stable(i)) then
               rbdn(i) = rbup(i)
               zl(i,k) = zl(i,k-1) - (t1(i,k)+t1(i,k-1))/2 *
     &              log(sl(i,k)/sl(i,k-1)) / gor
               thekv = t1(i,k)/slk(i,k)*(1.+fv*max(q1(i,k),qmin))
               spdk2 = max(rcl*(u1(i,k)**2+v1(i,k)**2),1.)
               rbup(i) = (thekv-the1v(i))*(g*zl(i,k)/the1v(i))/spdk2
               kpbl(i) = k
               stable(i) = rbup(i).gt.rbcr
            endif
         enddo
      enddo
c
      do i = 1,im
         k = kpbl(i)
         if(rbdn(i).ge.rbcr) then
            rbint = 0.
         elseif(rbup(i).le.rbcr) then
            rbint = 1.
         else
            rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
         endif
         hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
         if(hpbl(i).lt.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
      enddo
c
      do i = 1, im
         hol = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if(sfcflg(i)) then
            hol = min(hol,-zfmin)
         else
            hol = max(hol,zfmin)
         endif
c
         hol = hol*hpbl(i)/zl1(i)*sfcfrac
         if(sfcflg(i)) then
            phim(i) = (1.-aphi16*hol)**(-1./4.)
            phih(i) = (1.-aphi16*hol)**(-1./2.)
         else
            phim(i) = (1.+aphi5*hol)
            phih(i) = phim(i)
         endif
         wscale(i) = ustar(i)/phim(i)
         wscale(i) = min(wscale(i),ustar(i)*aphi16)
         wscale(i) = max(wscale(i),ustar(i)/aphi5)
      enddo
c
c     compute the surface variables for pbl height estimation
c     under unstable conditions
c
      do i = 1,im
         sflux  = heat(i) + evap(i)*fv*the1(i)
         if(sfcflg(i).and.sflux.gt.0.0) then
           hgamt(i) = min(cfac*heat(i)/wscale(i),gamcrt)
           hgamq(i) = min(cfac*evap(i)/wscale(i),gamcrq)
#ifdef ISOTOPE
           hgamqi1(i)=min(cfac*evapi1(i)/wscale(i),gamcrq)
           hgamqi2(i)=min(cfac*evapi2(i)/wscale(i),gamcrq)
           if(slmask(i).ne.1) then
              hgamq(i) = 0.
              hgamqi1(i) = 0.
              hgamqi2(i) = 0.
           endif
#else
           if(slmask(i).ne.1) hgamq(i) = 0.
#endif
           vpert = hgamt(i) + fv*the1(i)*hgamq(i)
           vpert = min(vpert,gamcrt)
           thermal(i) = thermal(i) + max(vpert,0.)
           hgamt(i) = max(hgamt(i),0.0)
           hgamq(i) = max(hgamq(i),0.0)
#ifdef ISOTOPE
           hgamqi1(i) = max(hgamqi1(i),0.0)
           hgamqi2(i) = max(hgamqi2(i),0.0)
#endif
         else
           pblflg(i) = .false.
         endif
      enddo
c
      do i = 1,im
         if(pblflg(i)) then
            kpbl(i) = 1
            hpbl(i) = zi(i,2)
         endif
      enddo
c
c     enhance the pbl height by considering the thermal
c
      do i = 1, im
         if(pblflg(i)) then
            stable(i) = .false.
            rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
         do i = 1, im
            if(.not.stable(i).and.pblflg(i)) then
               rbdn(i) = rbup(i)
               zl(i,k) = zl(i,k-1) - (t1(i,k)+t1(i,k-1))/2 *
     &              log(sl(i,k)/sl(i,k-1)) / gor
               thekv = t1(i,k)/slk(i,k)*(1.+fv*max(q1(i,k),qmin))
               spdk2 = max(rcl*(u1(i,k)**2+v1(i,k)**2),1.)
               rbup(i) = (thekv-thermal(i))*(g*zl(i,k)/the1v(i))/spdk2
               kpbl(i) = k
               stable(i) = rbup(i).gt.rbcr
            endif
         enddo
      enddo
c
      do i = 1,im
         if(pblflg(i)) then
            k = kpbl(i)
            if(rbdn(i).ge.rbcr) then
               rbint = 0.
            elseif(rbup(i).le.rbcr) then
               rbint = 1.
            else
               rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
            endif
            hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpbl(i).lt.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
            if(kpbl(i).le.1) pblflg(i) = .false.
         endif
      enddo
c
c     compute diffusion coefficients below pbl
c
      do k = 1, kmpbl
         do i = 1, im
            if(kpbl(i).gt.k) then
               prnum = (phih(i)/phim(i)+cfac*vk*.1)
               prnum = min(prnum,prmax)
               prnum = max(prnum,prmin)
               zfac = max((1.-(zi(i,k+1)-zl1(i))/
     1              (hpbl(i)-zl1(i))), zfmin)
               dku(i,k) = xkzo + wscale(i)*vk*zi(i,k+1)
     1              *zfac**pfac
               dkt(i,k) = dku(i,k)/prnum
               dku(i,k) = min(dku(i,k),dkmax)
               dku(i,k) = max(dku(i,k),dkmin)
               dkt(i,k) = min(dkt(i,k),dkmax)
               dkt(i,k) = max(dkt(i,k),dkmin)
            endif
         enddo
      enddo
c
c     compute diffusion coefficients over pbl (free atmosphere)
c
      do k = 1, km-1
         do i = 1, im
            if(k.ge.kpbl(i)) then
               ti   =0.5*(t1(i,k)+t1(i,k+1))
               rdz  =rdzt(i,k)/ti
               dw2  =rcl*((u1(i,k)-u1(i,k+1))**2+(v1(i,k)-v1(i,k+1))**2)
               shr2 =max(dw2,dw2min)*rdz**2
               tvd  =t1(i,k)*(1.+fv*max(q1(i,k),qmin))
               tvu  =t1(i,k+1)*(1.+fv*max(q1(i,k+1),qmin))
               bvf2 =g*(gocp+rdz*(tvu-tvd))/ti
               ri   =max(bvf2/shr2,rimin)
               zk   =vk*zi(i,k+1)
               rl2  =(zk*rlam/(rlam+zk))**2
               dk   =rl2*sqrt(shr2)
               if(ri.lt.0.) then ! unstable regime
                  sri = sqrt(-ri)
                  dku(i,k) = xkzo + dk*(1+8.*(-ri)/(1+1.746*sri))
                  dkt(i,k) = xkzo + dk*(1+8.*(-ri)/(1+1.286*sri))
               else             ! stable regime
                  dkt(i,k)  = xkzo + dk/(1+5.*ri)**2
                  prnum     = 1.0 + 2.1*ri
                  prnum     = min(prnum,prmax)
                  dku(i,k)  = (dkt(i,k)-xkzo)*prnum + xkzo
               endif
c
               dku(i,k) = min(dku(i,k),dkmax)
               dku(i,k) = max(dku(i,k),dkmin)
               dkt(i,k) = min(dkt(i,k),dkmax)
               dkt(i,k) = max(dkt(i,k),dkmin)
c
            endif
         enddo
      enddo
c
c     compute tridiagonal matrix elements for heat and moisture
c
      do i=1,im
         ad(i,1) = 1.
#ifndef VIC
         a1(i,1) = t1(i,1)-betat(i)*(t1(i,1)/slk(i,1)-tsea(i))
#ifdef HYDRO
         a2(i,1) = q1(i,1)-betaq(i)*(q1(i,1) - qss(i))
#else
         a2(i,1) = q1(i,1)-betaq(i)*(max(q1(i,1),qmin)-max(qss(i),qmin))
#ifdef ISOTOPE
         dq1=(max(q1(i,1),qmin)-max(qss(i),qmin))
#ifdef NOFRAC
         rsea(i,1)=1.D0
         rsea(i,2)=0.D0
#else
         rsea(i,1)=1.D0
         rsea(i,2)=1.D0
#endif
         if (dq1.le.0.D0) then !! evaporation
#ifdef NOFRAC
         if (xlon(i)/3.1416*180.gt.180 .and.
     $       xlon(i)/3.1416*180.lt.240 .and.
     $       xlat(i)/3.1416*180.gt.10  .and.
     $       xlat(i)/3.1416*180.lt.30 ) then
           revap(i,1)=1.D0 !! sea
           revap(i,2)=0.D0 !! land
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
           flg2(i)=2
         elseif (xlon(i)/3.1416*180.gt.220 .and.
     $       xlon(i)/3.1416*180.lt.250 .and.
     $       xlat(i)/3.1416*180.gt.30  .and.
     $       xlat(i)/3.1416*180.lt.50 ) then
           revap(i,1)=0.D0 !! sea
           revap(i,2)=1.D0 !! land
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
           flg2(i)=3
         else
           revap(i,1)=0.D0 !! sea
           revap(i,2)=0.D0 !! land
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
         endif
#else
         if (slmask(i).eq.0.D0) then
         if ((q1(i,1+2*km).gt.0.D0).and.
     $        (q1(i,1).gt.0.D0).and.
     $        (qssi(i,1).gt.0.D0).and.
     $        (qss(i).gt.0.D0).and.
     $        (qssi(i,2).gt.0.D0).and.
     $        (q1(i,1+3*km).gt.0.D0)) then            
           revap(i,1)=frkin(1,spd1(i))
           revap(i,2)=frkin(2,spd1(i))
           a2(i,1+2*km) = q1(i,1+2*km)-revap(i,1)*betaq(i)*
     $              (max(q1(i,1+2*km),qmin*q1(i,1+2*km)/q1(i,1))
     $              -max(qssi(i,1),qmin*qssi(i,1)/qss(i)))
           a2(i,1+3*km) = q1(i,1+3*km)-revap(i,2)*betaq(i)*
     $              (max(q1(i,1+3*km),qmin*q1(i,1+3*km)/q1(i,1))
     $              -max(qssi(i,2),qmin*qssi(i,2)/qss(i)))
           flg2(i)=1
         else
           revap(i,1)=1.D0/freq(1,tskin(i))*rsea(i,1)
           revap(i,2)=1.D0/freq(2,tskin(i))*rsea(i,2)
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
           flg2(i)=2
         endif
         elseif (slmask(i).ge.1.D0) then !! land and ice
           revap(i,1)=rresv1(i)
           revap(i,2)=rresv2(i)
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
           flg2(i)=3
         endif
#endif
         else !! dew
           if (q1(i,1).gt.0.D0) then
           revap(i,1)=freq(1,tskin(i))*q1(i,1+2*km)/q1(i,1)
           revap(i,2)=freq(2,tskin(i))*q1(i,1+3*km)/q1(i,1)
           a2(i,1+2*km)=q1(i,1+2*km)-
     $          revap(i,1)*betaq(i)*dq1
           a2(i,1+3*km)=q1(i,1+3*km)-
     $          revap(i,2)*betaq(i)*dq1
           flg2(i)=4
           else
           revap(i,1)=0.D0
           revap(i,2)=0.D0
           a2(i,1+2*km)=q1(i,1+2*km)
           a2(i,1+3*km)=q1(i,1+3*km)
           flg2(i)=5
           endif
         endif
#endif
#endif
#else
         a1(i,1) = t1(i,1)+beta2(i)*heat(i)
         a2(i,1) = q1(i,1)+beta2(i)*evap(i)
#endif
#ifdef CO2PROG
** co2 **
         a2(i,1+2*km)=q1(i,1+2*km)+dt*rdzt1(i)/t1(i,1)*co2fl(i)
#endif
      enddo
      if(ndiff.ge.2) then
#ifndef CO2PROG
#ifdef ISOTOPE
        do k = 2, 2
#else
        do k = 2, ndiff
#endif
#else
        do k = 2, 2
#endif
          is = (k-1) * km
          do i = 1, im
            a2(i,1+is) = q1(i,1+is)
          enddo
        enddo
      endif
c
      do k = 1,km-1
         do i = 1,im
            dtodsd = dt/del(i,k)
            dtodsu = dt/del(i,k+1)
            dsig   = sl(i,k)-sl(i,k+1)
            rdz = rdzt(i,k)*2./(t1(i,k)+t1(i,k+1))
            if(pblflg(i).and.k.lt.kpbl(i)) then
               dsdzt = dsig*dkt(i,k)*rdz*(gocp-hgamt(i)/hpbl(i))
               dsdzq = dsig*dkt(i,k)*rdz*(-hgamq(i)/hpbl(i))
               a2(i,k)   = a2(i,k)+dtodsd*dsdzq
               a2(i,k+1) = q1(i,k+1)-dtodsu*dsdzq
#ifdef ISOTOPE
               dsdzq = dsig*dkt(i,k)*rdz*(-hgamqi1(i)/hpbl(i))
               a2(i,k+2*km)   = a2(i,k+2*km)+dtodsd*dsdzq
               a2(i,k+1+2*km) = q1(i,k+1+2*km)-dtodsu*dsdzq
               dsdzq = dsig*dkt(i,k)*rdz*(-hgamqi2(i)/hpbl(i))
               a2(i,k+3*km)   = a2(i,k+3*km)+dtodsd*dsdzq
               a2(i,k+1+3*km) = q1(i,k+1+3*km)-dtodsu*dsdzq
#endif
            else
               dsdzt = dsig*dkt(i,k)*rdz*(gocp)
               a2(i,k+1) = q1(i,k+1)
#ifdef ISOTOPE
               a2(i,k+1+2*km) = q1(i,k+1+2*km)
               a2(i,k+1+3*km) = q1(i,k+1+3*km)
#endif
            endif
            dsdz2 = dsig*dkt(i,k)*rdz*rdz
            au(i,k)   = -dtodsd*dsdz2
            al(i,k)   = -dtodsu*dsdz2
            ad(i,k)   = ad(i,k)-au(i,k)
            ad(i,k+1) = 1.-al(i,k)
            a1(i,k)   = a1(i,k)+dtodsd*dsdzt
            a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
         enddo
      enddo
      if(ndiff.ge.2) then
#ifdef ISOTOPE
        do kk = 2, 2
#else
        do kk = 2, ndiff
#endif
          is = (kk-1) * km
          do k = 1, km - 1
            do i = 1, im
              a2(i,k+1+is) = q1(i,k+1+is)
            enddo
          enddo
        enddo
      endif
c
c     solve tridiagonal problem for heat and moisture
c
      call tridin(im,ILOTS,km,ndiff,al,ad,au,a1,a2,au,a1,a2)
c
c     recover tendencies of heat and moisture
c
      do  k = 1,km
         do i = 1,im
            ttend = (a1(i,k)-t1(i,k))*rdt
            qtend = (a2(i,k)-q1(i,k))*rdt
            tau(i,k) = tau(i,k)+ttend
            rtg(i,k) = rtg(i,k)+qtend
            dtsfc(i) = dtsfc(i)+cont*del(i,k)*pstar(i)*ttend
            dqsfc(i) = dqsfc(i)+conq*del(i,k)*pstar(i)*qtend
#ifdef ISOTOPE
            qtend = (a2(i,k+2*km)-q1(i,k+2*km))*rdt
            rtg(i,k+2*km) = rtg(i,k+2*km)+qtend
            dqsfci1(i) = dqsfci1(i)+conq*del(i,k)*pstar(i)*qtend
            qtend = (a2(i,k+3*km)-q1(i,k+3*km))*rdt
            rtg(i,k+3*km) = rtg(i,k+3*km)+qtend
            dqsfci2(i) = dqsfci2(i)+conq*del(i,k)*pstar(i)*qtend
#endif
#ifdef CO2PROG
** co2 **
            qtend = (a2(i,k+2*km)-q1(i,k+2*km))*rdt
            rtg(i,k+2*km) = rtg(i,k+2*km)+qtend
            dqsfci1(i) = dqsfci1(i)+conq*del(i,k)*pstar(i)*qtend
#endif
         enddo
      enddo
      if(ndiff.ge.2) then
#ifndef CO2PROG
#ifdef ISOTOPE
        do kk = 2, 2
#else
        do kk = 2, ndiff
#endif
#else
        do kk = 2, 2
#endif
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (a2(i,k+is)-q1(i,k+is))*rdt
              rtg(i,k+is) = rtg(i,k+is)+qtend
            enddo
          enddo
        enddo
      endif
c
c     compute tridiagonal matrix elements for momentum
c
      do i = 1,im
         ad(i,1) = 1.+betaw(i)
         a1(i,1) = u1(i,1)
         a2(i,1) = v1(i,1)
      enddo
c
      do k = 1,km-1
         do i=1,im
            dtodsd = dt/del(i,k)
            dtodsu = dt/del(i,k+1)
            dsig   = sl(i,k)-sl(i,k+1)
            rdz   = rdzt(i,k)*2./(t1(i,k)+t1(i,k+1))
            dsdz2 = dsig*dku(i,k)*rdz*rdz
            au(i,k)  = -dtodsd*dsdz2
            al(i,k)  = -dtodsu*dsdz2
            ad(i,k)  = ad(i,k)-au(i,k)
            ad(i,k+1)= 1.-al(i,k)
            a1(i,k+1)= u1(i,k+1)
            a2(i,k+1)= v1(i,k+1)
         enddo
      enddo
c
c     solve tridiagonal problem for momentum
c
      call tridi2(im,ILOTS,km,al,ad,au,a1,a2,au,a1,a2)
c
c     recover tendencies of momentum
c
      conwrc = conw*sqrt(rcl)
      do k = 1,km
         do i = 1,im
            utend = (a1(i,k)-u1(i,k))*rdt
            vtend = (a2(i,k)-v1(i,k))*rdt
            b(i,k)= b(i,k)+utend
            a(i,k)= a(i,k)+vtend
            dusfc(i) = dusfc(i)+conwrc*del(i,k)*pstar(i)*utend
            dvsfc(i) = dvsfc(i)+conwrc*del(i,k)*pstar(i)*vtend
         enddo
      enddo
c
      return
      end
