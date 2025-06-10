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
      subroutine cloud3(im,ix,kx,deltim,ps,
     1     t1,q1,q2,ncloud,sl,del,slk,rain,lat,dot,kdt,thour)
!-------------------------------------------------
!
!
!  This code is a simple ice microphyiscs scheme (WSM3) of the WRF
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
!           Jimy Dudhia (NCAR) and Shua-Hua Chen (UC Davis)
!           Summer 2003
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Dudhia (D89, 1989) J. Atmos. Sci.
!
#include <paramodel.h>
!
#undef NON
cc#define DEBUG
!
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
!     passing variables
      dimension del(im,kx),sl(im,kx),slk(im,kx),
     $          ps(im),kuo(im),dot(im,kx)
      dimension q1(ix,kx),q2(ix,kx,ncloud),t1(ix,kx),rain(im)
      dimension qrs(ix,kx)
!     local variables
!
!  local variables and arrays
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
      dimension psfc(ILOTS)
      dimension q(ILOTS,levs_),t(ILOTS,levs_)
      dimension qci(ILOTS,levs_)
      dimension rh(ILOTS,levs_),qs(ILOTS,levs_),p(ILOTS,levs_)
      dimension den(ILOTS,levs_),denfac(ILOTS,levs_)
      dimension rslope(ILOTS,levs_),rslope2(ILOTS,levs_),
     1          rslope3(ILOTS,levs_),rslopeb(ILOTS,levs_),
     1          pgen(ILOTS,levs_),paut(ILOTS,levs_),
     1          pacr(ILOTS,levs_),
     2          pisd(ILOTS,levs_),pres(ILOTS,levs_),
     1          pcon(ILOTS,levs_),
     2          fall(ILOTS,levs_),falk(ILOTS,levs_),
     3          mstep(ILOTS)
      dimension xl(ILOTS,levs_),cpm(ILOTS,levs_),
     1          work1(ILOTS,levs_),work2(ILOTS,levs_),
     2          xni(ILOTS,levs_)
#ifdef NON
      dimension delz(ILOTS,levs_)
#endif
      dimension numdt(ILOTS)
      logical   flgcld(ILOTS),qminus
!     physical parameters
!     parameter(g=#g,rd=#rd,rv=#rv,t0c=#t0c)
      parameter(g=g_,rd=rd_,rv=rv_,t0c=273.15)
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
!===============================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!     tvcal(x,y) = x+x*delta*max(y,qmin)
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-30 kg/kg.
!     otherwise use uniform distribution value (1.e38)
!
      lamdar(x,y)=(pidn0r/(x*y))**.25
      lamdas(x,y,z)=(pidn0s*z/(x*y))**.25
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
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = kts, kte
        do i = its, ite
          cpm(i,k) = cpmcal(q1(i,k))
          xl(i,k) = xlcal(t1(i,k))
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
      do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      enddo
!
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
!     initialize the large scale variables
!
      do k = kts, kte
        do i = its, ite
          p(i,k) = psfc(i) * sl(i,k)
          t(i,k) = t1(i,k)
          q(i,k) = q1(i,k)
          qci(i,k) = max(q2(i,k,1),0.)
          qrs(i,k) = max(q2(i,k,2),0.)
          den(i,k) = p(i,k)/(rd*t1(i,k))
          denfac(i,k) = sqrt(den0/den(i,k))
#ifdef NON
          delz(i,k) = psfc(i)*del(i,k)/g/den(i,k)
#endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          qs(i,k) = 1000.*fpvs(t(i,k))
          qs(i,k) = eps * qs(i,k) / (p(i,k)+epsm1*qs(i,k))
          qs(i,k) = max(qs(i,k),qmin)
          rh(i,k) = max(q(i,k) / qs(i,k),qmin)
#ifdef DEBUG
          qdt = dtcld * 1000.
          if(lat.eq.latd.and.i.eq.lond) write(iun,700) qci(i,k)*1000.,
     1      qrs(i,k)*1000.,qci(i,k)*1000.,qrs(i,k)*1000.
700       format(1x,'before computation', 4f10.4)
#endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = kts, kte
        do i = its, ite
          pres(i,k) = 0.
          paut(i,k) = 0.
          pacr(i,k) = 0.
          pgen(i,k) = 0.
          pisd(i,k) = 0.
          pcon(i,k) = 0.
          fall(i,k) = 0.
          falk(i,k) = 0.
#ifdef VCI
          fallc(i,k) = 0.
          falkc(i,k) = 0.
#endif
          xni(i,k) = 1.e3
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
!
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
#ifdef N0S_CONST
          n0sfac(i,k) = 1.
#else
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
#endif
          if(t(i,k).ge.t0c) then
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopermax
              rslopeb(i,k) = rsloperbmax
              rslope2(i,k) = rsloper2max
              rslope3(i,k) = rsloper3max
            else
              rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
              rslopeb(i,k) = rslope(i,k)**bvtr
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          else
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopesmax
              rslopeb(i,k) = rslopesbmax
              rslope2(i,k) = rslopes2max
              rslope3(i,k) = rslopes3max
            else
              rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
              rslopeb(i,k) = rslope(i,k)**bvts
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          endif
#ifdef FLETCHER
            xni(i,k) = max(xn0 * exp(betai*supcol),5.)
#else
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
            xni(i,k) = min(max(5.38e7*(den(i,k)                         
     1                *max(qci(i,k),qmin))**0.75,1.e3),1.e6)
#endif
        enddo
      enddo
!
#ifdef DEBUG
      if(lat.eq.latd) write(iun,605)
#endif
      mstepmax = 1
      numdt = 1
      do k = kte, kts, -1
        do i = its, ite
          if(t(i,k).lt.t0c) then
            pvt = pvts
          else
            pvt = pvtr
          endif
          work1(i,k) = pvt*rslopeb(i,k)*denfac(i,k)
#ifdef NON
          work2(i,k) = work1(i,k)/delz(i,k)
#else
          work2(i,k) = g*den(i,k)/psfc(i)*work1(i,k)/del(i,k)
#endif
          numdt(i) = max(nint(work2(i,k)*dtcld+.5),1)
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
            falk(i,k) = den(i,k)*qrs(i,k)*work2(i,k)/mstep(i)
            hold = falk(i,k)
            fall(i,k) = fall(i,k)+falk(i,k)
            holdrs = qrs(i,k)
            qrs(i,k) = max(qrs(i,k)-falk(i,k)*dtcld/den(i,k),0.)
#ifdef DEBUG
              if(lat.eq.latd.and.i.eq.lond)                              
     1          write(iun,606) k,p(i,k)/100.,                            
     2          t(i,k)-t0c,q(i,k)*1000.,rh(i,k)*100.,dot(i,k)*10000.,             
     3          work1(i,k),hold*qdt,falk(i,k)*qdt,holdrs*1000.,          
     4          qrs(i,k)*1000.,n
#endif
          endif
        enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if(n.le.mstep(i)) then
              falk(i,k) = den(i,k)*qrs(i,k)*work2(i,k)/mstep(i)
              hold = falk(i,k)
              fall(i,k) = fall(i,k)+falk(i,k)
              holdrs = qrs(i,k)
#ifdef NON
              qrs(i,k) = max(qrs(i,k)-(falk(i,k)
     1                  -falk(i,k+1)*delz(i,k+1)/delz(i,k))*dtcld
     2                  /den(i,k),0.)
#else
              qrs(i,k) = max(qrs(i,k)-(falk(i,k)
     1                  -falk(i,k+1)*del(i,k+1)/del(i,k))*dtcld
     2                  /den(i,k),0.)
#endif
#ifdef DEBUG
              if(lat.eq.latd.and.i.eq.lond)                              
     1          write(iun,606) k,p(i,k)/100.,                            
     2          t(i,k)-t0c,q(i,k)*1000.,rh(i,k)*100.,dot(i,k)*10000.,
     3          work1(i,k),hold*qdt,falk(i,k)*qdt,holdrs*1000.,          
     4          qrs(i,k)*1000.,n
#endif
            endif
          enddo
        enddo
      enddo
#ifdef DEBUG
605   format(1x,'  k','     p','    t','    q','   rh','     w',
     1      '    vt','  falk','  falt','  qrsi','  qrsf',' mstep')
606   format(1x,i3,f6.0,2f5.1,f5.0,f6.2,5f6.2,i5)
#endif
!
#ifdef VCI
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      mstepmax = 1
      do i = its, ite
        numdt(i) = 1
        mstep = 1
      enddo
      do k = kte, kts, -1
        do i = its, ite
          if(t(i,k).lt.t0c.and.qci(i,k).gt.0.) then
#ifdef HD1990
            work1c(i,k) = 3.29*(den(i,k)*qci(i,k))**0.16
#else
            xmi = den(i,k)*qci(i,k)/xni(i,k)
            diameter  = dicon * sqrt(xmi)
            work1c(i,k) = 1.49e4*diameter**1.31
#endif
          else
            work1c(i,k) = 0.
          endif
          if(qci(i,k).le.0.) then
            work2c(i,k) = 0.
          else
#ifdef NON
            work2c(i,k) = work1c(i,k)/delz(i,k)
#else
            work2c(i,k) = g*den(i,k)/psfc(i)*work1c(i,k)/del(i,k)
#endif
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
          if (n.le.mstep(i)) then
            falkc(i,k) = den(i,k)*qci(i,k)*work2c(i,k)/mstep(i)
            holdc = falkc(i,k)
            fallc(i,k) = fallc(i,k)+falkc(i,k)
            holdci = qci(i,k)
            qci(i,k) = max(qci(i,k)-falkc(i,k)*dtcld/den(i,k),0.)
          endif
        enddo
        do k = kte-1, kts, -1
          do i = its, ite
            if (n.le.mstep(i)) then
              falkc(i,k) = den(i,k)*qci(i,k)*work2c(i,k)/mstep(i)
              holdc = falkc(i,k)
              fallc(i,k) = fallc(i,k)+falkc(i,k)
              holdci = qci(i,k)
#ifdef NON
              qci(i,k) = max(qci(i,k)-(falkc(i,k)
     1                  -falkc(i,k+1)*delz(i,k+1)/delz(i,k))*dtcld
     2                  /den(i,k),0.)
#else
              qci(i,k) = max(qci(i,k)-(falkc(i,k)
     1                  -falkc(i,k+1)*del(i,k+1)/del(i,k))*dtcld
     2                  /den(i,k),0.)
#endif
            endif
          enddo
        enddo
      enddo
#endif
!
!----------------------------------------------------------------
!     compute the freezing/melting term. [D89 B16-B17]
!     freezing occurs one layer above the melting level
!
      do i = its, ite
        mstep(i) = 0
      enddo
      do k = kts, kte
!
        do i = its, ite
          if(t(i,k).ge.t0c) then
            mstep(i) = k
          endif
        enddo
      enddo
!
      do i = its, ite
        if(mstep(i).ne.0.and.dot(i,mstep(i)).lt.0.) then
          work1(i,1) = float(mstep(i) + 1)
          work1(i,2) = float(mstep(i))
        else
          work1(i,1) = float(mstep(i))
          work1(i,2) = float(mstep(i))
        endif
      enddo
!
      do i = its, ite
        k  = nint(work1(i,1))
        kk = nint(work1(i,2))
        if(k*kk.ge.1) then
          qrsci = qrs(i,k) + qci(i,k)
          if(qrsci.gt.0..or.fall(i,kk).gt.0.) then
#ifdef NON
            frzmlt = min(max(-w(i,k)*qrsci/delz(i,k),-qrsci/dtcld),
     1               qrsci/dtcld)
#else
            frzmlt = min(max(dot(i,k)*1000.*qrsci/del(i,k)/psfc(i),
     1               -qrsci/dtcld),qrsci/dtcld)
#endif
            snomlt = min(max(fall(i,kk)/den(i,kk),-qrs(i,k)/dtcld),    
     1               qrs(i,k)/dtcld)
            if(k.eq.kk) then
              t(i,k) = t(i,k) - xlf0/cpm(i,k)*(frzmlt+snomlt)*dtcld
            else
              t(i,k) = t(i,k) - xlf0/cpm(i,k)*frzmlt*dtcld
              t(i,kk) = t(i,kk) - xlf0/cpm(i,kk)*snomlt*dtcld
            endif

#ifdef DEBUG
            if(lat.eq.latd.and.i.eq.lond) write(iun,608) k,t(i,k)-t0c, 
     1        -dot(i,k)*10000,frzmlt*qdt,snomlt*qdt
#endif
          endif
        endif
      enddo
#ifdef DEBUG
 608  format(1x,'k = ',i3,' t = ',f5.1,' w = ',f6.2,' frz/mlt = ',f5.1, 
     1           ' snomlt = ',f5.1)
#endif
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
      do i = its, ite
        if(fall(i,1).gt.0.) then
#ifdef NON
          rain(i) = fall(i,1)*delz(i,1)/denr*dtcld
     1            + rain(i)
#else
          rain(i) = fall(i,1)*psfc(i)*del(i,1)/g/den(i,1)/denr*dtcld
     1            + rain(i)
#endif
        endif
      enddo
#ifdef DEBUG
      if(lat.eq.latd) write(iun,601) latd,lond,loop,rain(lond)
 601  format(1x,' ncloud3 lat lon loop : rain(mm) ',3i6,f20.2)
#endif
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
          if(t(i,k).ge.t0c) then
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopermax
              rslopeb(i,k) = rsloperbmax
              rslope2(i,k) = rsloper2max
              rslope3(i,k) = rsloper3max
            else
              rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
              rslopeb(i,k) = rslope(i,k)**bvtr
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          else
            if(qrs(i,k).le.qcrmin)then
              rslope(i,k) = rslopesmax
              rslopeb(i,k) = rslopesbmax
              rslope2(i,k) = rslopes2max
              rslope3(i,k) = rslopes3max
            else
              rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
              rslopeb(i,k) = rslope(i,k)**bvts
              rslope2(i,k) = rslope(i,k)*rslope(i,k)
              rslope3(i,k) = rslope2(i,k)*rslope(i,k)
            endif
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          if(t(i,k).ge.t0c) then
            work1(i,k) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k))
          else
            work1(i,k) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k))
          endif
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          supsat = max(q(i,k),qmin)-qs(i,k)
          satdt = supsat/dtcld
          if(t(i,k).ge.t0c) then
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!---------------------------------------------------------------
! paut1: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
!
            if(qci(i,k).gt.qc0) then
#ifdef KESSLER
              paut(i,k) = qck1*(qci(i,k)-qc0)
#else
              paut(i,k) = qck1*qci(i,k)**(7./3.)
#endif
              paut(i,k) = min(paut(i,k),qci(i,k)/dtcld)
            endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [D89 B15]
!        (C->R)
!---------------------------------------------------------------
            if(qrs(i,k).gt.qcrmin.and.qci(i,k).gt.qmin) then
                pacr(i,k) = min(pacrr*rslope3(i,k)*rslopeb(i,k) 
     1               *qci(i,k)*denfac(i,k),qci(i,k)/dtcld)
            endif
!---------------------------------------------------------------
! pres1: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
            if(qrs(i,k).gt.0.) then
                coeres = rslope2(i,k)*sqrt(rslope(i,k)*rslopeb(i,k))
                pres(i,k) = (rh(i,k)-1.)*(precr1*rslope2(i,k)            
     1                   +precr2*work2(i,k)*coeres)/work1(i,k)
              if(pres(i,k).lt.0.) then
                pres(i,k) = max(pres(i,k),-qrs(i,k)/dtcld)
                pres(i,k) = max(pres(i,k),satdt/2)
              else
                pres(i,k) = min(pres(i,k),satdt/2)
              endif
            endif
          else
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
            supcol = t0c-t(i,k)
            ifsat = 0
#ifdef FLETCHER
            xni(i,k) = max(xn0 * exp(betai*supcol),5.)
#else
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
            xni(i,k) = min(max(5.38e7*(den(i,k)                         
     1                *max(qci(i,k),qmin))**0.75,1.e3),1.e6)
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
            if(qrs(i,k).gt.qcrmin.and.qci(i,k).gt.qmin) then
#ifdef NEW_PSACI
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter  = min(dicon * sqrt(xmi),dimax)
#ifdef HD1990
              vt2i = 3.29*(den(i,k)*qci(i,k))**0.16
#else
              vt2i = 1.49e4*diameter**1.31
#endif
              vt2s = pvts*rslopeb(i,k)*denfac(i,k)
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k)+2.*diameter*rslope2(i,k)          
     1                +diameter**2*rslope(i,k)
              pacr(i,k) = min(pi*qci(i,k)*eacrs*n0s*n0sfac(i,k)          
     1                 *abs(vt2s-vt2i)*acrfac/4.,qci(i,k)/dtcld)
#else
              pacr(i,k) = min(pacrs*n0sfac(i,k)*eacrs*rslope3(i,k)       
     1               *rslopeb(i,k)*qci(i,k)*denfac(i,k),qci(i,k)/dtcld)
#endif
            endif
!-------------------------------------------------------------
! pisd: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
              xmi = den(i,k)*qci(i,k)/xni(i,k)
              diameter = dicon * sqrt(xmi)
              pisd(i,k) = 4.*diameter*xni(i,k)*(rh(i,k)-1.)              
     1                  /work1(i,k)
              if(pisd(i,k).lt.0.) then
                pisd(i,k) = max(pisd(i,k),satdt/2)
                pisd(i,k) = max(pisd(i,k),-qci(i,k)/dtcld)
              else
                pisd(i,k) = min(pisd(i,k),satdt/2)
              endif
              if(abs(pisd(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pres2: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k)*sqrt(rslope(i,k)*rslopeb(i,k))
              pres(i,k) = (rh(i,k)-1.)*n0sfac(i,k)*(precs1*rslope2(i,k)   
     1                  +precs2*work2(i,k)*coeres)/work1(i,k)
              if(pres(i,k).lt.0.) then
                pres(i,k) = max(pres(i,k),-qrs(i,k)/dtcld)
                pres(i,k) = max(pres(i,k),satdt/2)
              else
                pres(i,k) = min(pres(i,k),satdt/2)
              endif
              if(abs(pisd(i,k)+pres(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pgen: generation(nucleation) of ice from vapor [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
#ifdef NEW_ICE_GEN
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pgen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k),0.))
     1                   /dtcld)
#else
              pgen(i,k) = max(0.,(xm0*xni(i,k)/den(i,k)-max(qci(i,k)     
     1                   ,0.))/dtcld)
#endif
              pgen(i,k) = min(pgen(i,k),satdt)
            endif
!-------------------------------------------------------------
! paut2: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k).gt.0.) then
#ifdef NEW_ICE_AUTO
              qimax = roqimax/den(i,k)
              paut(i,k) = max(0.,(qci(i,k)-qimax)/dtcld)
#else
              paut(i,k) = max(0.,(qci(i,k)-xmmax*xni(i,k)
     1                   /den(i,k))/dtcld)
#endif
            endif
#ifdef DEBUG
      qdt = dtcld * 1000.
      if(lat.eq.latd.and.i.eq.lond) write(iun,801) k,xm0,xni(i,k),
     1   satdt,n0sfac(i,k),eacrs,rslope2(i,k),rslopeb(i,k),pacrs,qimax
801   format(1x,2i,9e11.3)
#endif
          endif
        enddo
      enddo
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = kts, kte
        do i = its, ite
          qciik = max(qmin,qci(i,k))
          delqci = (paut(i,k)+pacr(i,k)-pgen(i,k)-pisd(i,k))*dtcld
          if(delqci.ge.qciik) then
            facqci = qciik/delqci
            paut(i,k) = paut(i,k)*facqci
            pacr(i,k) = pacr(i,k)*facqci
            pgen(i,k) = pgen(i,k)*facqci
            pisd(i,k) = pisd(i,k)*facqci
          endif
          qik = max(qmin,q(i,k))
          delq = (pres(i,k)+pgen(i,k)+pisd(i,k))*dtcld
          if(delq.ge.qik) then
            facq = qik/delq
            pres(i,k) = pres(i,k)*facq
            pgen(i,k) = pgen(i,k)*facq
            pisd(i,k) = pisd(i,k)*facq
          endif
          work2(i,k) = -pres(i,k)-pgen(i,k)-pisd(i,k)
          q(i,k) = q(i,k)+work2(i,k)*dtcld
          qci(i,k) = max(qci(i,k)-(paut(i,k)+pacr(i,k)-pgen(i,k) 
     1             -pisd(i,k))*dtcld,0.)
          qrs(i,k) = max(qrs(i,k)+(paut(i,k)+pacr(i,k)               
     1             +pres(i,k))*dtcld,0.)
          if(t(i,k).lt.t0c) then
            t(i,k) = t(i,k)-xls*work2(i,k)/cpm(i,k)*dtcld
          else
            t(i,k) = t(i,k)-xl(i,k)*work2(i,k)/cpm(i,k)*dtcld
          endif
        enddo
      enddo
!
      do k = kts, kte
        do i = its, ite
          qs(i,k) = 1000.*fpvs0(t(i,k))
          qs(i,k) = eps * qs(i,k) / (p(i,k)+epsm1*qs(i,k))
          qs(i,k) = max(qs(i,k),qmin)
          denfac(i,k) = sqrt(den0/den(i,k))
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
          work1(i,k) = conden(t(i,k),q(i,k),qs(i,k),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k)+work1(i,k)
          pcon(i,k) = min(max(work1(i,k),0.),max(q(i,k),0.))/dtcld
          if(qci(i,k).gt.0..and.work1(i,k).lt.0.and.t(i,k).gt.t0c)   
     1      pcon(i,k) = max(work1(i,k),-qci(i,k))/dtcld
          q(i,k) = q(i,k)-pcon(i,k)*dtcld
          qci(i,k) = max(qci(i,k)+pcon(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcon(i,k)*xl(i,k)/cpm(i,k)*dtcld
#ifdef DEBUG
          if(lat.eq.latd.and.i.eq.lond) then
            pgens(i,k) = pgens(i,k)+pgen(i,k)
            pcons(i,k) = pcons(i,k)+pcon(i,k)
            pisds(i,k) = pisds(i,k)+pisd(i,k)
            pacrss(i,k) = pacrss(i,k)+pacr(i,k)
            press(i,k) = press(i,k)+pres(i,k)
            pauts(i,k) = pauts(i,k)+paut(i,k)
            write(iun,604) k,p(i,k)/100.,                                   
     1        t(i,k)-t0c,t(i,k)-t1(i,k),q(i,k)*1000.,                       
     2        (q(i,k)-q1(i,k))*1000.,rh(i,k)*100.,pgens(i,k)*qdt,           
     3        pcons(i,k)*qdt,pisds(i,k)*qdt,pauts(i,k)*qdt,pacrss(i,k)
     4        *qdt,press(i,k)*qdt,qci(i,k)*1000.,qrs(i,k)*1000.
          endif
#endif
        enddo
      enddo
#ifdef DEBUG
603   format(1x,'  k','     p',                                           
     1     '    t',' delt','    q',' delq','   rh',                       
     2     ' pgen',' pcon',' pisd',' paut',' pacr',' pres',               
     3     '  qci','  qrs')
604   format(1x,i3,f6.0,4f5.1,f5.0,8f5.2)
#endif
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = kts, kte
         do i = its, ite
            t1(i,k) = t(i,k)
            q1(i,k) = q(i,k)
            q2(i,k,1) = q2(i,k,1)+qci(i,k)-max(q2(i,k,1),qmin)
            q2(i,k,2) = q2(i,k,2)+qrs(i,k)-max(q2(i,k,2),qmin)
         enddo
      enddo
!
      enddo                  ! big loops
      return
      end
