#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define LEVS levs
#define LONF lonf
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONFP lonfp
#define RD rd
#define RV rv
#define G g
#define PI pi
!------------------------------------------------------------
      subroutine cloud6(im,ix,kx,deltim,ps,
     1     t1,q1,q2,ncloud,sl,del,slk,rain,lat,dot,kdt,thour)
!-------------------------------------------------
!
!  This code is a GRAUPEL phase ice microphyiscs scheme (WSM6) of the WRF
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. Related changes in ice-microphysics and description of
!  other microphysics are described in Hong et al. (2004).
!  all units are m.k.s. and source/sink terms are kgkg-1s-1.
!
! WRFSMMP cloud scheme
!
!  Coded by Song-You Hong and Jeong-Ock Lim (Yonsei Univ.)
!           Jimy Dudhia (NCAR) and Shu-Hua Chen (UC Davis)
!           Summer 2003
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Lim (2004) Master thesis, Yonsei Univ.
!             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
!             Rutledge, Hobbs (RH, 1983) J. Atmos. Sci.
!             Rutledge, Hobbs (RH, 1984) J. Atmos. Sci.
!
#include <paramodel.h>
!
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
!     passing variables
      dimension del(im,kx),sl(im,kx),slk(im,kx),
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
      parameter (ice=3)
      dimension psfc(ILOTS)
      dimension q(ILOTS,levs_),t(ILOTS,levs_)
      dimension qci(ILOTS,levs_,ice),qrs(ILOTS,levs_,ice)
      dimension rh(ILOTS,levs_,ice),qs(ILOTS,levs_,ice)
      dimension p(ILOTS,levs_)
      dimension den(ILOTS,levs_),denfac(ILOTS,levs_)
      dimension rslope(ILOTS,levs_,ice),rslope2(ILOTS,levs_,ice),
     1          rslopeb(ILOTS,levs_,ice),rslope3(ILOTS,levs_,ice),
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
      parameter(g=g_,rd=rd_,rv=rv_,t0c=273.15,pi=pi_)
      parameter(eps=rd/rv,delta=1./eps-1.,epsm1=eps-1.)
      parameter(iun=84)
      parameter(dtcldcr=120.)
      parameter(qminus=.false.)
      real lamdar,lamdas,lamdag
      real n0r,n0s,n0smax,n0g
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
      common /concld6/pidn0g,bvtg,precg1,precg2,pvtg,qs0,n0g,
     1                dens,pacrg,avtr,g6pbr,cl
!
      parameter( dicon = 11.9 )
      parameter( peaut =.55, xmyu =1.718e-5, r0=.8e-5)
      parameter( xncr =3.e8)
!
      real eacrs, xni0, roqi0
      real n0sfac(ILOTS,levs_)
      dimension falkc(ILOTS,levs_), work1c(ILOTS,levs_),
     1          work2c(ILOTS,levs_), fallc(ILOTS,levs_)
      real holdc, holdci
      dimension pracw(ILOTS,levs_), psacw(ILOTS,levs_),
     1          pgacw(ILOTS,levs_), pgacr(ILOTS,levs_),
     2          pgacs(ILOTS,levs_), psaci(ILOTS,levs_),
     3          pgml (ILOTS,levs_), praci(ILOTS,levs_),
     4          piacr(ILOTS,levs_), pracs(ILOTS,levs_),
     5          psacr(ILOTS,levs_), pgaci(ILOTS,levs_),
     6          pseml(ILOTS,levs_), pgeml(ILOTS,levs_)
      dimension pgev(ILOTS,levs_)
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
      lamdas(x,y,z)=(pidn0s*z/(x*y))**.25
      lamdag(x,y)=(pidn0g/(x*y))**.25
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
      qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu
      qc0  = 4./3.*pi*denr*r0**3*xncr/den0
!===============================================================
!     compute preliminary variables
!
      delt = deltim * 2.
!

      its = 1
      ite = im
      kts = 1
      kte = kx
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
          qrs(i,k,3) = max(q2(i,k,5),qmin)
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
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = 0.
          work1(i,k,2) = 0.
          work1(i,k,3) = 0.
          pres(i,k,1) = 0.
          pres(i,k,2) = 0.
          pres(i,k,3) = 0.
          paut(i,k,1) = 0.
          paut(i,k,2) = 0.
          paut(i,k,3) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pgen(i,k) = 0.
          pisd(i,k) = 0.
          pcon(i,k) = 0.
          psml(i,k) = 0.
          pgml(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psev(i,k) = 0.
          pgev(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          xni(i,k) = 1.e3
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
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
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
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
          xni(i,k) = min(max(5.38e7*(den(i,k)                           
     1              *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
        enddo
      enddo
!
      mstepmax = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          work1(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)/delz(i,k)
          work1(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)/delz(i,k)
          work1(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)/delz(i,k)
          numdt(i) = max(nint(max(work1(i,k,1),work1(i,k,2)
     1             ,work1(i,k,3)) 
     1              *dtcld+.5),1)
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
              falk(i,k,3) = den(i,k)*qrs(i,k,3)*work1(i,k,3)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              fall(i,k,3) = fall(i,k,3)+falk(i,k,3)
              qrs(i,k,1) = max(qrs(i,k,1)-falk(i,k,1)*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-falk(i,k,2)*dtcld/den(i,k),0.)
              qrs(i,k,3) = max(qrs(i,k,3)-falk(i,k,3)*dtcld/den(i,k),0.)
            endif
          enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              falk(i,k,1) = den(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
              falk(i,k,2) = den(i,k)*qrs(i,k,2)*work1(i,k,2)/mstep(i)
              falk(i,k,3) = den(i,k)*qrs(i,k,3)*work1(i,k,3)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              fall(i,k,3) = fall(i,k,3)+falk(i,k,3)
              qrs(i,k,1) = max(qrs(i,k,1)-(falk(i,k,1)-falk(i,k+1,1)    
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-(falk(i,k,2)-falk(i,k+1,2)    
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
              qrs(i,k,3) = max(qrs(i,k,3)-(falk(i,k,3)-falk(i,k+1,3)    
     1                    *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
        do k = kte, kts, -1
          do i = its, ite
            if(n.le.mstep(i).and.t(i,k).gt.t0c) then
!---------------------------------------------------------------
! psml: melting of snow [RH83 A25]
!       (T>T0: S->R)
!---------------------------------------------------------------
              xlf = xlf0
              work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
              if(qrs(i,k,2).gt.0.) then
                coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)
     1                 *rslopeb(i,k,2))
                psml(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2. 
     1                     *n0sfac(i,k)*(precs1*rslope2(i,k,2)          
     1                     +precs2*work2(i,k)*coeres)
                psml(i,k) = min(max(psml(i,k)*dtcld/mstep(i),           
     1                      -qrs(i,k,2)/mstep(i)),0.)
                qrs(i,k,2) = qrs(i,k,2) + psml(i,k)
                qrs(i,k,1) = qrs(i,k,1) - psml(i,k)
                t(i,k) = t(i,k) + xlf/cpm(i,k)*psml(i,k)
              endif
!---------------------------------------------------------------
! pgml: melting of graupel [LFO 47]
!       (T>T0: G->R)
!---------------------------------------------------------------
              if(qrs(i,k,3).gt.0.) then
                coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)
     1                 *rslopeb(i,k,3))
                pgml(i,k) = xka(t(i,k),den(i,k))/xlf                    
     1                     *(t0c-t(i,k))*(precg1*rslope2(i,k,3)         
     1                     +precg2*work2(i,k)*coeres)
                pgml(i,k) = min(max(pgml(i,k)*dtcld/mstep(i),           
     1                      -qrs(i,k,3)/mstep(i)),0.)
                qrs(i,k,3) = qrs(i,k,3) + pgml(i,k)
                qrs(i,k,1) = qrs(i,k,1) - pgml(i,k)
                t(i,k) = t(i,k) + xlf/cpm(i,k)*pgml(i,k)
              endif
            endif
          enddo
        enddo
      enddo
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
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = min(dicon * sqrt(xmi),dimax)
            work1c(i,k) = 1.49e4*diameter**1.31
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
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
      do i = its, ite
        fallsum = fall(i,kts,1)+fall(i,kts,2)+fall(i,kts,3)
        if(fallsum.gt.0.) then
          rain(i) = fallsum*delz(i,kts)/denr*dtcld + rain(i)
        endif
      enddo
!
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
          if(supcol.gt.0..and.qci(i,k,1).gt.qmin) then
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
     1            *(exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2             
     1            *rslope(i,k,1)*dtcld,qrs(i,k,1))
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
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
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k)
     1                 ,qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
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
            paut(i,k,1) = qck1*qci(i,k,1)**(7./3.)
            paut(i,k,1) = min(paut(i,k,1),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)        
     1                  *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
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
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
          xni(i,k) = min(max(5.38e7*(den(i,k)                           
     1                 *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          eacrs = exp(0.05*(-supcol))
!
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          if(supcol.gt.0.and.qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              praci(i,k) = pacrr*rslope3(i,k,1)*rslopeb(i,k,1)          
     1                    *qci(i,k,2)*denfac(i,k)
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
!-------------------------------------------------------------
! piacr: Accretion of rain by cloud ice [LFO 26]
!        (T<T0: R->S or R->G)
!-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)     
     1                    *g6pbr*rslope3(i,k,1)*rslope3(i,k,1)          
     1                    *rslopeb(i,k,1)/24./den(i,k)
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by snow [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin) then
              psaci(i,k) = pacrs*n0sfac(i,k)*eacrs*rslope3(i,k,2)       
     1                    *rslopeb(i,k,2)*qci(i,k,2)*denfac(i,k)
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! pgaci: Accretion of cloud ice by graupel [LFO 41]
!        (T<T0: I->G)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              pgaci(i,k) = pacrg*egi*rslope3(i,k,3)*rslopeb(i,k,3)      
     1                   *qci(i,k,2)*denfac(i,k)
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [LFO 24]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)           
     1                  *rslopeb(i,k,2)*qci(i,k,1)*denfac(i,k)          
     1                  ,qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacw: Accretion of cloud water by graupel [LFO 40]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)        
     1                  *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pracs: Accretion of snow by rain [LFO 27]
!         (T<T0: S->G)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)   
     1                +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)  
     1                +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2s)     
     1                    *(dens/den(i,k))*acrfac
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! psacr: Accretion of rain by snow [LFO 28]
!         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
!-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)     
     1              +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)    
     1              +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2s-vt2r)       
     1                  *(denr/den(i,k))*acrfac
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacr: Accretion of rain by graupel [LFO 42]
!         (T<T0: R->G) (T>=T0: enhance melting of graupel)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)     
     1              +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)    
     1              +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2g-vt2r)*(denr/den(i,k))   
     1                  *acrfac
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacs: Accretion of snow by graupel [LFO 29]
!        (S->G)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,2).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,3)     
     1              +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,3)    
     1              +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,3)
            if(supcol.gt.0) then
              egs = exp(-0.09*supcol)
            else
              egs = 1.
            endif
            pgacs(i,k) = pi**2*egs*n0s*n0sfac(i,k)*n0g*abs(vt2g-vt2s)   
     1                  *(dens/den(i,k))*acrfac
            pgacs(i,k) = min(pgacs(i,k),qrs(i,k,2)/dtcld)
          endif
          if(supcol.le.0) then
            xlf = xlf0
!-------------------------------------------------------------
! pseml: Enhanced melting of snow by accretion of water
!        (T>=T0: S->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.)                                        
     1        pseml(i,k) = min(max(cl*supcol*(psacw(i,k)+psacr(i,k))  
     1                    /xlf,-qrs(i,k,2)/dtcld),0.)
!-------------------------------------------------------------
! pgeml: Enhanced melting of graupel by accretion of water [RH84 A21-A22]
!        (T>=T0: G->R)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0.)                                        
     1        pgeml(i,k) = min(max(cl*supcol*(pgacw(i,k)+pgacr(i,k))  
     1                    /xlf,-qrs(i,k,3)/dtcld),0.)
          endif
          if(supcol.gt.0) then
!-------------------------------------------------------------
! pisd: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
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
!-------------------------------------------------------------
! pres2: deposition/sublimation rate of snow [HDC 14]
!        (T<T0: V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              pres(i,k,2) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1          
     1                     *rslope2(i,k,2)+precs2*work2(i,k)            
     1                     *coeres)/work1(i,k,2)
              supice = satdt-pres(i,k,1)-pisd(i,k)
              if(pres(i,k,2).lt.0.) then
                pres(i,k,2) = max(pres(i,k,2),-qrs(i,k,2)/dtcld)
                pres(i,k,2) = max(max(pres(i,k,2),satdt/2),supice)
              else
                pres(i,k,2) = min(min(pres(i,k,2),satdt/2),supice)
              endif
              if(abs(pres(i,k,1)+pisd(i,k)+pres(i,k,2)).ge.abs(satdt))  
     1          ifsat = 1
            endif
!-------------------------------------------------------------
! pres3: deposition/sublimation rate of graupel [LFO 46]
!        (T<T0: V->G or G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pres(i,k,3) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)       
     1                        +precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-pres(i,k,1)-pisd(i,k)-pres(i,k,2)
              if(pres(i,k,3).lt.0.) then
                pres(i,k,3) = max(pres(i,k,3),-qrs(i,k,3)/dtcld)
                pres(i,k,3) = max(max(pres(i,k,3),satdt/2),supice)
              else
                pres(i,k,3) = min(min(pres(i,k,3),satdt/2),supice)
              endif
              if(abs(pres(i,k,1)+pisd(i,k)+pres(i,k,2)+pres(i,k,3)).ge. 
     1          abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pgen: generation(nucleation) of ice from vapor [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-pres(i,k,1)-pisd(i,k)-pres(i,k,2)
     1               -pres(i,k,3)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pgen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))    
     1                   /dtcld)
              pgen(i,k) = min(min(pgen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! paut2: conversion(aggregation) of ice to snow [HDC 12]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              paut(i,k,2) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif
!
!-------------------------------------------------------------
! paut3: conversion(aggregation) of snow to graupel [LFO 37]
!        (T<T0: S->G)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              paut(i,k,3) = min(max(0.,alpha2*(qrs(i,k,2)-qs0))         
     1                     ,qrs(i,k,2)/dtcld)
            endif
          endif
!
!-------------------------------------------------------------
! psev: Evaporation of melting snow [RH83 A27]
!       (T>=T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psev(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1            
     1                     *rslope2(i,k,2)+precs2*work2(i,k)            
     1                     *coeres)/work1(i,k,1)
              psev(i,k) = min(max(psev(i,k),-qrs(i,k,2)/dtcld),0.)
            endif
!-------------------------------------------------------------
! pgev: Evaporation of melting graupel [RH84 A19]
!       (T>=T0: G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgev(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)         
     1                   +precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgev(i,k) = min(max(pgev(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
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
!
          delta2=0.
          delta3=0.
          if(qrs(i,k,1).lt.1.e-4.and.qrs(i,k,2).lt.1.e-4) delta2=1.
          if(qrs(i,k,1).lt.1.e-4) delta3=1.
          if(t(i,k).le.t0c) then
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source = (paut(i,k,1)+pracw(i,k)+psacw(i,k)+pgacw(i,k))
     1             *dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
            endif
!
!     cloud ice
!
            value = max(qmin,qci(i,k,2))
            source = (paut(i,k,2)-pgen(i,k)-pisd(i,k)+praci(i,k)        
     1              +psaci(i,k)+pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,2) = paut(i,k,2)*factor
              pgen(i,k) = pgen(i,k)*factor
              pisd(i,k) = pisd(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif
!
!     rain
!
            value = max(qmin,qrs(i,k,1))
            source = (-paut(i,k,1)-pres(i,k,1)-pracw(i,k)+piacr(i,k)    
     1              +psacr(i,k)+pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pres(i,k,1) = pres(i,k,1)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif
!
!     snow
!
            value = max(qmin,qrs(i,k,2))
            source = -(pres(i,k,2)+paut(i,k,2)-paut(i,k,3)              
     1               +piacr(i,k)*delta3+praci(i,k)*delta3               
     1               -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2          
     1               +psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              pres(i,k,2) = pres(i,k,2)*factor
              paut(i,k,2) = paut(i,k,2)*factor
              paut(i,k,3) = paut(i,k,3)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
!
!     graupel
!
            value = max(qmin,qrs(i,k,3))
            source = -(pres(i,k,3)+paut(i,k,3)+psacw(i,k)               
     1               +piacr(i,k)*(1.-delta3)+praci(i,k)*(1.-delta3)     
     1               +psacr(i,k)*(1.-delta2)+pracs(i,k)*(1.-delta2)     
     1               +pgaci(i,k)+pgacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pres(i,k,3) = pres(i,k,3)*factor
              paut(i,k,3) = paut(i,k,3)*factor
              psacw(i,k) = psacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
!
            work2(i,k)=-(pres(i,k,1)+pres(i,k,2)+pres(i,k,3)+pgen(i,k)   
     1                 +pisd(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(paut(i,k,1)+pracw(i,k)          
     1                     +psacw(i,k)+pgacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(paut(i,k,1)+pracw(i,k)          
     1                     +pres(i,k,1)-piacr(i,k)-pgacr(i,k)            
     1                     -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(paut(i,k,2)+praci(i,k)          
     1                     +psaci(i,k)+pgaci(i,k)-pgen(i,k)-pisd(i,k))   
     1                     *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(pres(i,k,2)+paut(i,k,2)         
     1                     -paut(i,k,3)+piacr(i,k)*delta3                
     1                     +praci(i,k)*delta3+psaci(i,k)-pgacs(i,k)      
     1                     -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2)    
     1                     *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pres(i,k,3)+paut(i,k,3)         
     1                     +psacw(i,k)+piacr(i,k)*(1.-delta3)            
     1                     +praci(i,k)*(1.-delta3)+psacr(i,k)
     1                     *(1.-delta2)
     1                     +pracs(i,k)*(1.-delta2)+pgaci(i,k)+pgacw(i,k) 
     1                     +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(pres(i,k,2)+pres(i,k,3)+pisd(i,k)+pgen(i,k)) 
     1                -xl(i,k)*pres(i,k,1)-xlf*(piacr(i,k)+psacw(i,k)    
     1                +pgacw(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source=(paut(i,k,1)+pracw(i,k)+psacw(i,k)+pgacw(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
            endif
!
!     rain
!
            value = max(qmin,qrs(i,k,1))
            source = (-psacw(i,k)-paut(i,k,1)+pseml(i,k)+pgeml(i,k)     
     1               -pracw(i,k)-pgacw(i,k)-pres(i,k,1))*dtcld
            if (source.gt.value) then
              factor = value/source
              paut(i,k,1) = paut(i,k,1)*factor
              pres(i,k,1) = pres(i,k,1)*factor
              pgacw(i,k) = pgacw(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
!
!     snow
!
            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psev(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psev(i,k) = psev(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif
!
!     graupel
!
            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgev(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgev(i,k) = pgev(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
            work2(i,k)=-(pres(i,k,1)+psev(i,k)+pgev(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(paut(i,k,1)+pracw(i,k)         
     1              +psacw(i,k)+pgacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(paut(i,k,1)+pracw(i,k)         
     1              +pres(i,k,1)+psacw(i,k)+pgacw(i,k)-pseml(i,k)       
     1              -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psev(i,k)-pgacs(i,k)           
     1              +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgev(i,k)           
     1              +pgeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(pres(i,k,1)+psev(i,k)+pgev(i,k))        
     1                -xlf*(pseml(i,k)+pgeml(i,k))
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
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
          qs(i,k,2) = 1000.*fpvs(t(i,k))
          qs(i,k,2) = eps * qs(i,k,2) / (p(i,k)+epsm1*qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcon: condensational/evaporational rate of cloud water [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
      do k = kts, kte
        do i = its, ite
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k)
     1                  ,cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcon(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)
     1              /dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                   
     1      pcon(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcon(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcon(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcon(i,k)*xl(i,k)/cpm(i,k)*dtcld
        enddo
      enddo
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
          q2(i,k,5) = q2(i,k,5)+qrs(i,k,3)-max(q2(i,k,5),qmin)
        enddo
      enddo
!
      enddo                  ! big loops

      return
      end
