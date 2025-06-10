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
#define RD rd
#define RV rv
!------------------------------------------------------------
      subroutine cloud5(im,ix,kx,deltim,ps,
     1     t1,q1,q2,ncloud,sl,del,slk,rain,lat,dot,kdt,thour)
!-------------------------------------------------
!
!  This code is a MIXED phase ice microphyiscs scheme (WSM5) of the WRF
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. Related changes in ice-microphysics and description of
!  other microphysics are described in Hong et al. (2004).
!  all units are m.k.s. and source/sink terms are kgkg-1s-1.
!
! WRFSMMP cloud scheme
!
!  Coded by Song-You Hong, Jeong-Ock Lim (Yonsei Univ.)
!           Jimy Dudhia (NCAR) and Shu-Hua Chen (UC Davis)
!           Summer 2003
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Rutledge, Hobbs (RH, 1983) J. Atmos. Sci.
!
#include <paramodel.h>
!
#undef DEBUG
!
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
!     passing variables
      dimension del(ix,kx),sl(ix,kx),slk(ix,kx),
     $          ps(im),kuo(im),dot(im,kx)
      dimension q1(ix,kx),q2(ix,kx,ncloud),t1(ix,kx),rain(im)
!     local variables
!
!  local variables and arrays
!
!
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
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
      parameter (ice=2)
      dimension psfc(ILOTS)
      dimension q(ILOTS,levs_),t(ILOTS,levs_)
      dimension qci(ILOTS,levs_,ice),qrs(ILOTS,levs_,ice)
      dimension rh(ILOTS,levs_,ice),qs(ILOTS,levs_,ice)
      dimension p(ILOTS,levs_)
      dimension den(ILOTS,levs_),denfac(ILOTS,levs_)
      dimension rslope(ILOTS,levs_,ice),rslope2(ILOTS,levs_,ice),
     1          rslope3(ILOTS,levs_,ice),
     1          rslopeb(ILOTS,levs_,ice),
     1          pgen(ILOTS,levs_),paut(ILOTS,levs_,ice),
     1          pacr(ILOTS,levs_,3),
     2          pisd(ILOTS,levs_),pres(ILOTS,levs_,ice),
     1          pcon(ILOTS,levs_),psml(ILOTS,levs_),
     1          psev(ILOTS,levs_),
     2          fall(ILOTS,levs_,ice),falk(ILOTS,levs_,ice),
     3          mstep(ILOTS)
      dimension xl(ILOTS,levs_),cpm(ILOTS,levs_),
     1          work1(ILOTS,levs_,ice),work2(ILOTS,levs_),
     2          xni(ILOTS,levs_)
      dimension delz(ILOTS,levs_)
      dimension numdt(ILOTS)
      logical   flgcld(ILOTS),qminus
!     physical parameters
!     parameter(g=#g,rd=#rd,rv=#rv,t0c=#t0c)
      parameter(g=g_,rd=rd_,rv=rv_,t0c=273.15)
      parameter(pi=pi_)
      parameter(eps=rd/rv,delta=1./eps-1.,epsm1=eps-1.)
      parameter(iun=84)
      parameter(dtcldcr=120.)
      parameter(qminus=.false.)
      real lamdar,lamdas
      real n0r,n0s,n0smax
      common /concld35/den0,denr,cpd,cpv,xlv0,xlv1,xlf0,xls,
     1             sigtop,pidn0r,pidn0s,alpha,n0r,n0s,
     2             bvtr,bvtr2,bvtr3,pvtr,pacrr,precr1,precr2,
     3  betai,xn0,xm0,xmmax,bvts,bvts2,bvts3,pvts,pacrs,precs1,precs2,
     4  t40c,pacrc,pfrz1,pfrz2,qmin,qcrmin,
     5  n0smax,dimax,roqimax,
     6       rslopermax,rslopesmax,rslopegmax,       
     7       rsloperbmax,rslopesbmax,rslopegbmax,   
     8       rsloper2max,rslopes2max,rslopeg2max,  
     9       rsloper3max,rslopes3max,rslopeg3max
!
#ifdef NEW_ICE_GEN
      parameter( dicon = 11.9 )
#else
      parameter( dicon = 16.3 )
#endif
#ifdef KESSLER
      parameter( qck1 = 1.e-3,qc0 = 0.5e-3)
#else
      parameter( peaut =.55, xmyu =1.718e-5, r0=.8e-5)
#endif
      parameter( xncr =3.e8)
!
#ifdef NEW_MPHYS
      real eacrs, xni0, roqi0
#endif
      real n0sfac(ILOTS,levs_)
#ifdef VCI
      dimension falkc(ILOTS,levs_), work1c(ILOTS,levs_),
     1          work2c(ILOTS,levs_), fallc(ILOTS,levs_)
      real holdc, holdci
#endif
#ifdef DEBUG
      DIMENSION pgens(ILOTS,levs_),
     1          pauts(ILOTS,levs_), pacrss(ILOTS,levs_),
     2          pisds(ILOTS,levs_), press(ILOTS,levs_),
     3          pcons(ILOTS,levs_)
#endif
!
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!     tvcal(x,y) = x+x*ep1*max(y,qmin)
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
!
      lamdar(x,y)=(pidn0r/(x*y))**.25
#ifdef DEN_SNOW_OPT
      lamdas(x,y,z)=(pi*xcon*n0s*z/(x*y))**xlampow
#else
      lamdas(x,y,z)=(pidn0s*z/(x*y))**.25
#endif
!
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!
      diffus(x,y) = 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = (viscos(b,c)/diffus(b,a))**(.3333333)       
     1       /viscos(b,c)**(.5)*(den0/c)**0.25
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!
#ifndef KESSLER
      qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu
      qc0  = 4./3.*pi_*denr*r0**3*xncr/den0
#endif
!===============================================================
!     compute preliminary variables
!
      delt = deltim * 2.
#ifdef DEBUG
!=================================================================
!     set iprt = 0 for no unit fort.84 output
!
      iprt = 1
      if(iprt.eq.1) then
        lond = (im+1)/2 + 1
        latd = 42
        lond = 53
      else
        latd = jts
        lond = 1
      endif
#endif
!
      its = 1
      ite = im
      kts = 1
      kte = 1
!
      do k = 1, kx
        do i = its, ite
        if((sl(i,k).gt.sigtop).and.(k.gt.kte-1)) kte = k + 1
        enddo
      enddo
      kte = min(kte,kx)
!
      do i = its, ite
         psfc(i) = ps(i) * 1000.
         rain(i) = 0.
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = kts, kte
        do i = its, ite
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
        enddo
      enddo
#ifdef DEBUG
      if(lat.eq.latd) then
        i = lond
        do k = kts, kte
          press(i,k) = 0.
          pauts(i,k) = 0.
          pacrss(i,k)= 0.
          pgens(i,k) = 0.
          pisds(i,k) = 0.
          pcons(i,k) = 0.
          t1(i,k) = t(i,k)
          q1(i,k) = q(i,k)
        enddo
      endif
#endif
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo
!
!     remove negative qci qrs. this option will be taken out
!     when semi-lagrangian numeric scheme is introduced
!
      if(.not.qminus) then
      do n = 1,ncloud
      do k = kts, kte
        do i = its, ite
          q2(i,k,n) = max(q2(i,k,n),0.)
        enddo
      enddo
      enddo
      endif
!
!
!     initialize the large scale variables
!
      do k = kts, kte
        do i = its, ite
          p(i,k) = psfc(i) * sl(i,k)
          t(i,k) = t1(i,k)
          q(i,k) = q1(i,k)
          qci(i,k,1) = max(q2(i,k,1),qmin)
          qci(i,k,2) = max(q2(i,k,2),qmin)
          qrs(i,k,1) = max(q2(i,k,3),qmin)
          qrs(i,k,2) = max(q2(i,k,4),qmin)
          den(i,k) = p(i,k)/(rd*t1(i,k))
          denfac(i,k) = sqrt(den0/den(i,k))
          delz(i,k) = psfc(i)*del(i,k)/g/den(i,k)
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          qs(i,k,1) = 1000.*fpvs0(t(i,k))
          qs(i,k,1) = eps * qs(i,k,1) / (p(i,k)+epsm1*qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
          qs(i,k,2) = 1000.*fpvs(t(i,k))
          qs(i,k,2) = eps * qs(i,k,2) / (p(i,k)+epsm1*qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
#ifdef DEBUG
          qdt = delt * 1000.
          if(lat.eq.latd.and.i.eq.lond) write(iun,700) qci(i,k,1)*1000., 
     1      qrs(i,k,1)*1000.,qci(i,k,2)*1000.,qrs(i,k,2)*1000.
700       format(1x,'before computation', 4f10.4)
#endif
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
        do i = its, ite
          pres(i,k,1) = 0.
          pres(i,k,2) = 0.
          paut(i,k,1) = 0.
          paut(i,k,2) = 0.
          pacr(i,k,1) = 0.
          pacr(i,k,2) = 0.
          pacr(i,k,3) = 0.
          pgen(i,k) = 0.
          pisd(i,k) = 0.
          pcon(i,k) = 0.
          psml(i,k) = 0.
          psev(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
#ifdef VCI
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
#endif
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
#ifdef N0S_CONST
          n0sfac(i,k) = 1.
#else
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
#endif
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
#ifdef FLETCHER
          xni(i,k) = max(xn0 * exp(betai*supcol),5.)
#else
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
          xni(i,k) = min(max(5.38e7*(den(i,k)                           
     1              *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
#endif
        enddo
      enddo
!
      mstepmax = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          work1(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)/delz(i,k)
          work1(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)/delz(i,k)
          numdt(i) = max(nint(max(work1(i,k,1),
     1               work1(i,k,2))*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = its, ite
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = kte
        do i = its, ite
          if(n.le.mstep(i)) then
              falk(i,k,1) = den(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
              falk(i,k,2) = den(i,k)*qrs(i,k,2)*work1(i,k,2)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              qrs(i,k,1) = max(qrs(i,k,1)-falk(i,k,1)*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-falk(i,k,2)*dtcld/den(i,k),0.)
            endif
          enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              falk(i,k,1) = den(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
              falk(i,k,2) = den(i,k)*qrs(i,k,2)*work1(i,k,2)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              qrs(i,k,1) = max(qrs(i,k,1)-(falk(i,k,1)-falk(i,k+1,1)    
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-(falk(i,k,2)-falk(i,k+1,2)    
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
        do k = kte, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              if(t(i,k).gt.t0c.and.qrs(i,k,2).gt.0.) then
!----------------------------------------------------------------
! psml: melting of snow [RH83 A25]
!       (T>T0: S->R)
!----------------------------------------------------------------
                xlf = xlf0
                work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
                coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)
     1                  *rslopeb(i,k,2))
                psml(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2. 
     1                      *n0sfac(i,k)*(precs1*rslope2(i,k,2)+precs2  
     2                      *work2(i,k)*coeres)
                psml(i,k) = min(max(psml(i,k)*dtcld/mstep(i),           
     1                      -qrs(i,k,2)/mstep(i)),0.)
                qrs(i,k,2) = qrs(i,k,2) + psml(i,k)
                qrs(i,k,1) = qrs(i,k,1) - psml(i,k)
                t(i,k) = t(i,k) + xlf/cpm(i,k)*psml(i,k)
              endif
            endif
          enddo
        enddo
      enddo
#ifdef VCI
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      mstepmax = 1
      mstep = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          if(qci(i,k,2).le.0.) then
            work2c(i,k) = 0.
          else
#ifdef HD1990
            work1c(i,k) = 3.29*(den(i,k)*qci(i,k,2))**0.16
#else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = min(dicon * sqrt(xmi),dimax)
            work1c(i,k) = 1.49e4*diameter**1.31
#endif
            work2c(i,k) = work1c(i,k)/delz(i,k)
          endif
          numdt(i) = max(nint(work2c(i,k)*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = its, ite
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = kte
        do i = its, ite
          if(n.le.mstep(i)) then
            falkc(i,k) = den(i,k)*qci(i,k,2)*work2c(i,k)/mstep(i)
            holdc = falkc(i,k)
            fallc(i,k) = fallc(i,k)+falkc(i,k)
            holdci = qci(i,k,2)
            qci(i,k,2) = max(qci(i,k,2)-falkc(i,k)*dtcld/den(i,k),0.)
          endif
        enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              falkc(i,k) = den(i,k)*qci(i,k,2)*work2c(i,k)/mstep(i)
              holdc = falkc(i,k)
              fallc(i,k) = fallc(i,k)+falkc(i,k)
              holdci = qci(i,k,2)
              qci(i,k,2) = max(qci(i,k,2)-(falkc(i,k)-falkc(i,k+1)      
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
      enddo
#endif
!
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
      do i = its, ite
        fallsum = fall(i,1,1)+fall(i,1,2)
        if(fallsum.gt.0.) then
          rain(i) = fallsum*delz(i,1)/denr*dtcld + rain(i)
        endif
      enddo
!
#ifdef DEBUG
      if(lat.eq.latd) write(iun,601) latd,lond,loop,rain(lond)
 601  format(1x,' ncloud5 lat lon loop : rain(mm) ',3i6,f20.2)
#endif
!---------------------------------------------------------------
! piml: instantaneous melting of cloud ice [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = 0.
          endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c
!        (T<-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = 0.
          endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qci(i,k,1).gt.0.) then
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                  
     1         *den(i,k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif
!---------------------------------------------------------------
! pfrz: freezing of rain water [LFO 45]
!        (T<T0, R->S)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
            pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)             
     1            *(exp(pfrz2*supcol)-1.)*rslope(i,k,1)**7*dtcld,       
     2            qrs(i,k,1))
            qrs(i,k,2) = qrs(i,k,2) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
!----------------------------------------------------------------
!     rsloper: reverse of the slope parameter of the rain(m)
!     xka:    thermal conductivity of air(jm-1s-1k-1)
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
      do k = kts, kte
        do i = its, ite
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k)
     1                  ,qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!
      do k = kts, kte
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld
!---------------------------------------------------------------
! paut1: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          if(qci(i,k,1).gt.qc0) then
#ifdef KESSLER
            paut(i,k,1) = qck1*(qci(i,k,1)-qc0)
#else
            paut(i,k,1) = qck1*qci(i,k,1)**(7./3.)
            paut(i,k,1) = min(paut(i,k,1),qci(i,k,1)/dtcld)
#endif
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pacr(i,k,1) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)       
     1                   *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pres1: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.0.) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            pres(i,k,1) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)         
     1                   +precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(pres(i,k,1).lt.0.) then
              pres(i,k,1) = max(pres(i,k,1),-qrs(i,k,1)/dtcld)
              pres(i,k,1) = max(pres(i,k,1),satdt/2)
            else
              pres(i,k,1) = min(pres(i,k,1),satdt/2)
            endif
          endif
        enddo
      enddo
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0
#ifdef FLETCHER
          xni(i,k) = max(xn0 * exp(betai*supcol),5.)
#else
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
          xni(i,k) = min(max(5.38e7*(den(i,k)                           
     1                 *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
#endif
#ifdef NEW_ICE_ACR
#ifdef NEW_PSACI
          eacrs = exp(0.07*(-supcol))
#else
          eacrs = exp(0.05*(-supcol))
#endif
#else
          eacrs = 0.1
#endif
!
          if(supcol.gt.0) then
            if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,2).gt.qmin) then
#ifdef NEW_PSACI
              xmi = den(i,k)*qci(i,k,2)/xni(i,k)
              diameter  = min(dicon * sqrt(xmi),dimax)
#ifdef HD1990
              vt2i = 3.29*(den(i,k)*qci(i,k,2))**0.16
#else
              vt2i = 1.49e4*diameter**1.31
#endif
              vt2s = pvts*rslopeb(i,k,2)*denfac(i,k)
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)     
     1               +diameter**2*rslope(i,k,2)
              pacr(i,k,2) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)         
     1                    *abs(vt2s-vt2i)*acrfac/4.
#else
              pacr(i,k,2) = min(pacrs*n0sfac(i,k)*eacrs*rslope3(i,k,2)  
     1                    *rslopeb(i,k,2)*qci(i,k,2)*denfac(i,k)       
     1                    ,qci(i,k,2)/dtcld)
#endif
            endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
              pacr(i,k,3) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)        
     1                     *rslopeb(i,k,2)*qci(i,k,1)*denfac(i,k)       
     2                     ,qci(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! pisd: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              xmi = den(i,k)*qci(i,k,2)/xni(i,k)
              diameter = dicon * sqrt(xmi)
              pisd(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)
     1                  /work1(i,k,2)
              supice = satdt-pres(i,k,1)
              if(pisd(i,k).lt.0.) then
                pisd(i,k) = max(max(pisd(i,k),satdt/2),supice)
                pisd(i,k) = max(pisd(i,k),-qci(i,k,2)/dtcld)
              else
                pisd(i,k) = min(min(pisd(i,k),satdt/2),supice)
              endif
              if(abs(pres(i,k,1)+pisd(i,k)).ge.abs(satdt)) ifsat = 1
            endif
          endif
!-------------------------------------------------------------
! pres2: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
            coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
            pres(i,k,2) = (rh(i,k,2)-1.)*n0sfac(i,k)                    
     1                   *(precs1*rslope2(i,k,2)+precs2                 
     2                   *work2(i,k)*coeres)/work1(i,k,2)
            supice = satdt-pres(i,k,1)-pisd(i,k)
            if(pres(i,k,2).lt.0.) then
              pres(i,k,2) = max(pres(i,k,2),-qrs(i,k,2)/dtcld)
              pres(i,k,2) = max(max(pres(i,k,2),satdt/2),supice)
            else
              pres(i,k,2) = min(min(pres(i,k,2),satdt/2),supice)
            endif
            if(abs(pres(i,k,1)+pisd(i,k)+pres(i,k,2)).ge.abs(satdt))    
     1        ifsat = 1
          endif
!-------------------------------------------------------------
! pgen: generation(nucleation) of ice from vapor [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
          if(supcol.gt.0) then
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-pres(i,k,1)-pisd(i,k)-pres(i,k,2)
#ifdef NEW_ICE_GEN
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pgen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))    
     1                   /dtcld)
#else
              pgen(i,k) = max(0.,(xm0*xni(i,k)/den(i,k)-max(qci(i,k,2)  
     1                   ,0.))/dtcld)
#endif
              pgen(i,k) = min(min(pgen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! paut2: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
#ifdef NEW_ICE_AUTO
              qimax = roqimax/den(i,k)
              paut(i,k,2) = max(0.,(qci(i,k,2)-qimax)/dtcld)
#else
              paut(i,k,2) = max(0.,(qci(i,k,2)-xmmax*xni(i,k)/den(i,k))/dtcld)
#endif
            endif
          endif
!-------------------------------------------------------------
! psev: Evaporation of melting snow [RH83 A27]
!       (T>T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.)                    
     1        psev(i,k) = pres(i,k,2)*work1(i,k,2)/work1(i,k,1)
              psev(i,k) = min(max(psev(i,k),-qrs(i,k,2)/dtcld),0.)
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = kts, kte
        do i = its, ite
          if(t(i,k).le.t0c) then
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source = (paut(i,k,1)+pacr(i,k,1)+pacr(i,k,3))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pacr(i,k,1) = pacr(i,k,1)*factor
              pacr(i,k,3) = pacr(i,k,3)*factor
            endif
!
!     cloud ice
!
            value = max(qmin,qci(i,k,2))
            source = (paut(i,k,2)+pacr(i,k,2)-pgen(i,k)-pisd(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,2) = paut(i,k,2)*factor
              pacr(i,k,2) = pacr(i,k,2)*factor
              pgen(i,k) = pgen(i,k)*factor
              pisd(i,k) = pisd(i,k)*factor
            endif
!
            work2(i,k)=-(pres(i,k,1)+pres(i,k,2)+pgen(i,k)+pisd(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(paut(i,k,1)+pacr(i,k,1)         
     1                  +pacr(i,k,3))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(paut(i,k,1)+pacr(i,k,1)         
     1                  +pres(i,k,1))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(paut(i,k,2)+pacr(i,k,2)         
     1                  -pgen(i,k)-pisd(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(pres(i,k,2)+paut(i,k,2)         
     1                  +pacr(i,k,2)+pacr(i,k,3))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(pres(i,k,2)+pisd(i,k)+pgen(i,k))             
     1                -xl(i,k)*pres(i,k,1)-xlf*pacr(i,k,3)
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source=(paut(i,k,1)+pacr(i,k,1)+pacr(i,k,3))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pacr(i,k,1) = pacr(i,k,1)*factor
              pacr(i,k,3) = pacr(i,k,3)*factor
            endif
!
!     snow
!
            value = max(qcrmin,qrs(i,k,2))
            source=(-psev(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psev(i,k) = psev(i,k)*factor
            endif
            work2(i,k)=-(pres(i,k,1)+psev(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(paut(i,k,1)+pacr(i,k,1)        
     1                  +pacr(i,k,3))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(paut(i,k,1)+pacr(i,k,1)        
     1                  +pres(i,k,1) +pacr(i,k,3))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+psev(i,k)*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(pres(i,k,1)+psev(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          qs(i,k,1) = 1000.*fpvs0(t(i,k))
          qs(i,k,1) = eps * qs(i,k,1) / (p(i,k)+epsm1*qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          qs(i,k,2) = 1000.*fpvs(t(i,k))
          qs(i,k,2) = eps * qs(i,k,2) / (p(i,k)+epsm1*qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcon: condensational/evaporational rate of cloud water [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
#ifdef DEBUG
      if(lat.eq.latd) write(iun,603)
#endif
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k)
     1                 ,cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcon(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)
     1              /dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                   
     1      pcon(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcon(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcon(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcon(i,k)*xl(i,k)/cpm(i,k)*dtcld
#ifdef DEBUG
          if(lat.eq.latd.and.i.eq.lond) then
            pgens(i,k) = pgens(i,k)+pgen(i,k)
            pcons(i,k) = pcons(i,k)+pcon(i,k)
            pisds(i,k) = pisds(i,k)+pisd(i,k)
            pacrss(i,k) = pacrss(i,k)+pacr(i,k,1)+pacr(i,k,2)+pacr(i,k,3)
            press(i,k) = press(i,k)+pres(i,k,1)+pres(i,k,2)
            pauts(i,k) = pauts(i,k)+paut(i,k,1)+paut(i,k,2)
            write(iun,604) k,p(i,k)/100.,                                   
     1        t(i,k)-t0c,t(i,k)-t1(i,k),q(i,k)*1000.,                       
     2        (q(i,k)-q1(i,k))*1000.,rh(i,k,2)*100.,pgens(i,k)*qdt,         
     3        pcons(i,k)*qdt,pisds(i,k)*qdt,pauts(i,k)*qdt,pacrss(i,k)*qdt, 
     4        press(i,k)*qdt,qci(i,k,1)*1000.,qrs(i,k,1)*1000.,             
     5        qci(i,k,2)*1000.,qrs(i,k,2)*1000.
          endif
#endif
        enddo
      enddo
#ifdef DEBUG
603   format(1x,'  k','     p',                                             
     1     '    t',' delt','    q',' delq','   rh',                         
     2     ' pgen',' pcon',' pisd',' paut',' pacr',' pres',                 
     3     '   qc','   qr','   qi','   qs')
604   format(1x,i3,f6.0,4f5.1,f5.0,10f5.2)
#endif
!
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = kts, kte
        do i = its, ite
          t1(i,k) = t(i,k)
          q1(i,k) = q(i,k)
          q2(i,k,1) = q2(i,k,1)+qci(i,k,1)-max(q2(i,k,1),qmin)
          q2(i,k,2) = q2(i,k,2)+qci(i,k,2)-max(q2(i,k,2),qmin)
          q2(i,k,3) = q2(i,k,3)+qrs(i,k,1)-max(q2(i,k,3),qmin)
          q2(i,k,4) = q2(i,k,4)+qrs(i,k,2)-max(q2(i,k,4),qmin)
        enddo
      enddo
!
      enddo                  ! big loops
      return
      end
