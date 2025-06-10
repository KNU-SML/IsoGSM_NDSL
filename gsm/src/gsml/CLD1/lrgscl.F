#define CP cp
#define G g
#define HVAP hvap
#define RD rd
#define RV rv
#include <define.h>
      subroutine lrgscl(ims2,imx2,kmx,dt,ps,t,q,sl,del,slk,
#ifdef ISOTOPE
     +                  ntrac,rni,trni,trn,
#endif
#ifdef SAS_DIAG
     +                  dlt,dlq,dlh,
#endif
#ifdef CONDENS
     +                  qcilrg,qrslrg,
#endif
#ifndef RASV2
     +                  rn,lat)
#else
     +                  rn,lat,clw)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    lrgscl      calculate grid-scale precipitation
c   prgmmr: hua-lu pan       org: w/nmc23    date: 94-04-15
c
c abstract: calculates grid-scale condensation for one leap-frog
c   timestep, produces rain, and adjusts temperature and specific
c   humidity by wet-bulb process.  evaporation of part or all of the
c   rain may occur as it traverses unsaturated layers on the way down.
c   for conditionally unstable layers, a convective adjustment procedure
c   is applied to adjust to a uniform theta-e.
c
c program history log:
c   94-04-15  hua-lu pan
c   07-04-25  isotope
c
c usage:    call lrgscl(im,im2,km,dt,ps,t,q,sl,del,slk,rain,lat)
c
c   input argument list:
c     im       - integer number of points
c     im2      - real first dimension of t and q
c     km       - integer number of levels
c     dt       - real time step in seconds
c     ps       - real (im) surface pressure in kilopascals (cb)
c     t        - real (im2,km) current temperature in kelvin
c     q        - real (im2,km) current specific humidity in kg/kg
c     sl       - real (km) sigma values
c     del      - real (km) sigma layer thickness
c     slk      - real (km) sigma values ** kappa
c     lat      - integer latitude number
c
c   output argument list:
c     q        - real (im2,km) adjusted specific humidity in kg/kg
c     t        - real (im2,km) adjusted temperature in kelvin
c     rn       - real (im) large-scale rain in meters
c
c subprograms called:
c   fpvs     - function to compute saturation vapor pressure
c   fpkap    - function to compute p raised to the factor kappa
c   fthe     - function to compute theta-e
c   ftma     - function to compute temperature and moisture along a
c              moist adiabat
c
c remarks: the precipitation reaching the ground should be halved
c   before it is used to increment geshem, the running total.  this
c   prevents double-counting of geshem, which is incremented every
c   half leap-frog timestep.
c          the evaporation rate for falling precip is calculated
c   according to a method devised by e.kessler, in which a mean drop
c   surface area is obtained from the rainwater content.
c          ice is not considered.
c          functions fpvs,fpkap,fthe,ftma are inlined by fpp.
c
c attributes:
c   language: fortran 77
c   machine:  cray
c
c$$$
#include <paramodel.h>
#ifdef ISOTOPE
      dimension del(imx2,kmx),sl(imx2,kmx),ps(imx2),
     &          q(imx2,kmx*ntrac),t(imx2,kmx),
     $          rn(imx2),rni(imx2,ntrac-2),
     $          trn(imx2),trni(imx2,ntrac-2)
      integer   ntrac,kk,flg
      real      ff,alpha(2),tmpi
      real      freq            !! funcion
      real      f1,f2,eqf,alp0,alp20
#ifdef NOFRAC
      parameter (f1=0.99D0,f2=0.60D0,eqf=0.0D0) !! f1: max rh, f2: min rh
#else
      parameter (f1=0.99D0,f2=0.60D0,eqf=0.95D0) !! f1: max rh, f2: min rh
#endif
      real      eqlee !! lee's equilibrium calculation
      parameter (eqlee=0) !! 1:on, 0:off
      real      dif18o,difhdo   !! ratio of diffusion coef.
#ifdef NOFRAC
      parameter (dif18o=1.0D0,difhdo=1.0D0)
#else
      parameter (dif18o=1.02849D0,difhdo=1.02512D0)
#endif
      real      nn              !! degree of freedom
#ifdef NOFRAC
      parameter (nn=0.0)
#else
      parameter (nn=0.58)
#endif
      real      beta,gamma,rirat,virat,rnio,rh
#else
      dimension del(imx2,kmx),sl(imx2,kmx),ps(imx2),
     &          q(imx2,kmx),t(imx2,kmx),rn(imx2)
#endif
c
#ifdef RASV2
      dimension clw(imx2,kmx)
#endif
#ifdef SAS_DIAG
      dimension dlt(imx2,kmx), dlq(imx2,kmx), dlh(imx2,kmx)
#endif
#ifdef CONDENS
      dimension qcilrg(imx2,kmx), qrslrg(imx2,kmx)
#endif
c
c  physical parameters
c
      parameter(g=g_,rd=rd_,rv=rv_,
     &          cp=cp_,hvap=hvap_)
      parameter(elocp=hvap/cp,el2orc=hvap*hvap/(rv*cp),
     &          eps=rd/rv,epsm1=rd/rv-1.)
#ifndef ICE
      parameter(nx=7501)
      dimension tbpvs0(nx)
      common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
c-----------------------------------------------------------------------
#ifdef SAS_DIAG
      do k=kmx,1,-1
        do i=1,ims2
          dlt(i,k) = 0.0
          dlq(i,k) = 0.0
          dlh(i,k) = 0.0
        enddo
      enddo
#endif
c-----------------------------------------------------------------------
c
c  condense rain to wetbulb temperature if supersaturated
c  or evaporate rain using kessler parameterization.
c
      do i=1,ims2
        rn(i)=0.
#ifdef ISOTOPE
        trn(i)=0.
        do k=1,2
          rni(i,k)=0.
          trni(i,k)=0.
        enddo
#endif
      enddo
      do k=kmx,1,-1
        do i=1,ims2
#ifdef CONDENS
          qcilrg(i,k) = 0.0
          qrslrg(i,k) = 0.0
#endif
          dpovg=del(i,k)/g*ps(i)
#ifdef ICE
          es=fpvs(t(i,k))
#else
          xj1=min(max(c1xpvs0+c2xpvs0*t(i,k),1.),float(nx))
          jx1=min(xj1,nx-1.)
          fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
          es=fpvs01
c         es=fpvs0(t(i,k))
#endif
          qs=eps*es/(sl(i,k)*ps(i)+epsm1*es)
          qcond=(q(i,k)-qs)/(1.+el2orc*qs/t(i,k)**2)
#ifdef RASV2
          rn(i)=rn(i)+clw(i,k)
#endif
#ifdef ISOTOPE
          flg=0
          ff=0.D0
          do kk=3,ntrac
            if (kk.eq.3) then
              alpha(1)=freq(1,t(i,k))
            elseif (kk.eq.4) then
              alpha(2)=freq(2,t(i,k))
            endif
*** Super saturation Kinetic fractionation ***
* alpha_eff=alpha_eq*alpha_kin
* alpha_kin=S/(alpha_eq*D_iso/D*(S-1)+1)
* S=1-0.003*T [C] (Jouzel and Merlivat 1984)
            if (t(i,k).le.253.15) then
              if (kk.eq.3) then
                alpha(kk-2)=alpha(kk-2)*
     $                (1-0.003*(t(i,k)-273.15))/
     $                (alpha(kk-2)*dif18o*(-0.003*(t(i,k)-273.15))+1)
              elseif (kk.eq.4) then
                alpha(kk-2)=alpha(kk-2)*
     $                (1-0.003*(t(i,k)-273.15))/
     $                (alpha(kk-2)*difhdo*(-0.003*(t(i,k)-273.15))+1)
              endif
*** linear interpolation between -20 to 0 C ***
            elseif (t(i,k).le.273.15) then
              if (kk.eq.3) then
                alp0=freq(1,273.15)
                alp20=freq(1,253.15)*
     $               (1-0.003*(-20.))/
     $               (freq(1,253.15)*dif18o*(-0.003*(-20.))+1)
                alpha(kk-2)=
     $                (alp0*(t(i,k)-253.15)
     $                +alp20*(273.15-t(i,k)))/20.
              elseif (kk.eq.4) then
                alp0=freq(2,273.15)
                alp20=freq(2,253.15)*
     $               (1-0.003*(-20.))/
     $               (freq(2,253.15)*difhdo*(-0.003*(-20.))+1)
                alpha(kk-2)=
     $                (alp0*(t(i,k)-253.15)
     $                +alp20*(273.15-t(i,k)))/20.
              endif
            endif
          enddo
c
#endif
          if(qcond.gt.0.) then  !! condensation
            q(i,k)=q(i,k)-qcond
            t(i,k)=t(i,k)+qcond*elocp
            rn(i)=rn(i)+qcond*dpovg
#ifdef ISOTOPE
            trn(i)=trn(i)+qcond*dpovg
            flg=7
#endif
#ifdef SAS_DIAG
            dlq(i,k) = dlq(i,k) - qcond
            dlt(i,k) = dlt(i,k) + qcond*elocp
#endif
#ifdef CONDENS
            qcilrg(i,k) = qcond
#endif
#ifdef ISOTOPE
            if (q(i,k).ge.0.D0) then !! q > q' >= 0.
            flg=1
c
            if (q(i,k)+qcond.gt.0.D0) then
              ff=q(i,k)/(q(i,k)+qcond)
            else
              print*,'SomeWrgLrg1',i,k,q(i,k)+qcond
              ff=1.D0
            endif
c
            do kk=3,ntrac
*** open system ***
              rni(i,kk-2)=rni(i,kk-2)+
     $             q(i,k+(kk-1)*kmx)*(1.D0-ff**alpha(kk-2))*dpovg
              trni(i,kk-2)=trni(i,kk-2)+
     $             q(i,k+(kk-1)*kmx)*(1.D0-ff**alpha(kk-2))*dpovg
              q(i,k+(kk-1)*kmx)=q(i,k+(kk-1)*kmx)*ff**alpha(kk-2)
*** closed system ***
c              rni(i,kk-2)=rni(i,kk-2)+
c     $             q(i,k+(kk-1)*kmx)*(1-ff)*
c     $             (ff*(1.-1./alpha(kk-2))+1.)*dpovg
c              q(i,k+(kk-1)*kmx)=
c     $             q(i,k+(kk-1)*kmx)*   ff *
c     $             (ff*(1.-1./alpha(kk-2))+1./alpha(kk-2))
            enddo            
c
            elseif (q(i,k).lt.-qcond) then !! -q > -q' > 0.
            flg=2
c
            alpha(1)=1.D0
            alpha(2)=1.D0
            if (rn(i).gt.0.D0) then
              ff=(rn(i)-qcond*dpovg)/rn(i)
            else
              print*,'SomeWrgLrg2',i,k,rn(i)
              ff=1.D0
            endif
c
            do kk=3,ntrac
              q(i,k+(kk-1)*kmx)=q(i,k+(kk-1)*kmx)+
     $             rni(i,kk-2)*(1.D0-ff**(-1./alpha(kk-2)))/dpovg
              trni(i,kk-2)=trni(i,kk-2)-
     $             rni(i,kk-2)*(1.D0-ff**(-1./alpha(kk-2)))
              rni(i,kk-2)=rni(i,kk-2)*ff**(-1./alpha(kk-2))
            enddo
c
            else !! q>=0. q'<=0
            flg=3
c
            alpha(1)=1.D0
            alpha(2)=1.D0
            if (rn(i).gt.0.D0) then
              ff=(rn(i)+q(i,j)*dpovg)/rn(i)
            else
              print*,'SomeWrgLrg3',i,k,rn(i)
              ff=1.D0
            endif
c
            do kk=3,ntrac
              tmpi=q(i,k+(kk-1)*kmx)
              tmpi2=rni(i,kk-2)
              q(i,k+(kk-1)*kmx)=
     $             (rni(i,kk-2)+tmpi*dpovg)*
     $             (1.D0-ff**(-1./alpha(kk-2)))/dpovg
              rni(i,kk-2)=(rni(i,kk-2)+tmpi*dpovg)*
     $             ff**(-1./alpha(kk-2))
              trni(i,kk-2)=trni(i,kk-2)+rni(i,kk-2)-tmpi2
            enddo
c
            endif
c
#endif
          elseif(rn(i).gt.0.D0) then
            qevap=-qcond*(1.-exp(-0.32*sqrt(2.*dt*rn(i))))
            rnevap=min(qevap*dpovg,rn(i))
            q(i,k)=q(i,k)+rnevap/dpovg
            t(i,k)=t(i,k)-rnevap/dpovg*elocp
            rn(i)=rn(i)-rnevap
#ifdef ISOTOPE
            flg=8
#endif
#ifdef SAS_DIAG
            dlq(i,k) = dlq(i,k) + rnevap/dpovg
            dlt(i,k) = dlt(i,k) - rnevap/dpovg*elocp
#endif
#ifdef CONDENS
            qcilrg(i,k) = qcilrg(i,k) - rnevap
#endif
#ifdef ISOTOPE
*** raindrop evaporation ***
            if (rn(i).ge.0.D0) then !! qL > qL' >= 0.
            flg=4
c
            if (rn(i).eq.0.D0) then !! all evaporation
            do kk=3,ntrac
              q(i,k+(kk-1)*kmx)=q(i,k+(kk-1)*kmx)+
     $             rni(i,kk-2)/dpovg
              rni(i,kk-2)=0.D0
            enddo
            else
c
*** liquid equilibrium ***
            if (rn(i)+rnevap.gt.0.D0) then
              ff=rn(i)/(rn(i)+rnevap)
            else
              print*,'SomeWrgLrg4',i,k,rn(i)+rnevap
              ff=1.D0
            endif
c
*** Stewart 1975 ***
*   rni/rn=(rni/ro-gamma*qi/q)*ff**beta+gamma*qi/q
*   gamma=alpha_eq*RH/(alpha_eq*(D_iso/D)**n*(RH-1)+1)
*   beta=1/(alpha_eq*(D_iso/D)**n*(1-RH))-1
*   RH: relative humidity
            rh=min(max(f2,1.4*q(i,k)/qs),f1)
c            rh=0.95
            do kk=3,ntrac
c
            if (t(i,k).ge.273.15) then
            if (kk.eq.3) then
               gamma=alpha(kk-2)*rh/
     $              (alpha(kk-2)*dif18o**nn*(rh-1)+1.D0)
               beta=1.D0/(alpha(kk-2)*dif18o**nn*(1-rh))-1.D0
            elseif (kk.eq.4) then
               gamma=alpha(kk-2)*rh/
     $              (alpha(kk-2)*difhdo**nn*(rh-1)+1)
               beta=1.D0/(alpha(kk-2)*difhdo**nn*(1-rh))-1.D0
            endif
            else                   !! ice -> vap
               gamma=1.D0
               beta=0.D0
            endif
*** lee's equilibrium model ***
            if (eqlee.eq.1) then
              dz=rd*t(i,k)*del(i,k)/sl(i,k)/g
              call eqm_deg(kk-1,t(i,k),rn(i)*1000./dt,
     $             ps(i)*sl(i,k)*1000.,dz,eqf1)
            else
*** constant equilibrium rate ***
              eqf1=eqf
            endif
c
            rirat=rni(i,kk-2)/(rn(i)+rnevap)
            if (q(i,k)-rnevap/dpovg.gt.0.e0) then
               virat=q(i,k+(kk-1)*kmx)/(q(i,k)-rnevap/dpovg)
            else
               virat=0.
            endif
            rnio=rni(i,kk-2)
            rni(i,kk-2)=
     $           (rirat*(eqf1*ff**beta+1.D0-eqf1)
     $           +eqf1*(1.D0-ff**beta)*gamma*virat)
     $           *rn(i)
            q(i,k+(kk-1)*kmx)=
     $           q(i,k+(kk-1)*kmx)+(rnio-rni(i,kk-2))/dpovg
c
            enddo
c
            endif
c
            elseif (rn(i).lt.-rnevap) then !! -qL > -qL' > 0.
            flg=5
c
            alpha(1)=1.D0
            alpha(2)=1.D0
            if (q(i,k).gt.0.D0) then
              ff=(q(i,k)-rnevap/dpovg)/q(i,k)
            else
              print*,'SomeWrgLrg5',i,k,q(i,k)
              ff=1.D0
            endif
c
            do kk=3,ntrac
              rni(i,kk-2)=rni(i,kk-2)+
     $             q(i,k+(kk-1)*kmx)*(1.D0-ff**(-alpha(kk-2)))*dpovg
              q(i,k+(kk-1)*kmx)=q(i,k+(kk-1)*kmx)*ff**(-alpha(kk-2))
            enddo
c
            else !! qL >= 0. qL' <= 0.
            flg=6
c
            alpha(1)=1.D0
            alpha(2)=1.D0
            if (q(i,k).gt.0.D0) then
              ff=(q(i,k)+rn(i)/dpovg)/q(i,k)
            else
              print*,'SomeWrgLrg6',i,k,q(i,k)
              ff=1.D0
            endif
c
            do kk=3,ntrac
              tmpi=rni(i,kk-2)
              rni(i,kk-2)=
     $             (q(i,k+(kk-1)*kmx)+tmpi/dpovg)*
     $             (1.D0-ff**(-alpha(kk-2)))*dpovg
              q(i,k+(kk-1)*kmx)=
     $             (q(i,k+(kk-1)*kmx)+tmpi/dpovg)*ff**(-alpha(kk-2))
            enddo
c
            endif
#endif
          endif
c
#ifdef SAS_DIAG
          dlh(i,k) = cp*dlt(i,k) + hvap*dlq(i,k)
#endif
#ifdef CONDENS
          qcilrg(i,k) = max(qcilrg(i,k),0.0)
#endif
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end
