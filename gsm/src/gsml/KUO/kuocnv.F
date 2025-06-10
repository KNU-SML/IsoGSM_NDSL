#define CP cp
#define G g
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define RD rd
#define RV rv
#define T0C t0c
      subroutine kuocnv(ims,im,km,jcap,dt,del,sl,slk,ps,qn,
     &                  q,t,rn,kbot,ktop,kuo)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpkap,fpvs,fpvs0)
#endif
cinclude "dbkuocnv"
#include <paramodel.h>
      dimension del(km),sl(km),slk(km),ps(im),qn(ims,km),
     &      q(ims,km),t(ims,km),rn(im),kbot(im),ktop(im),kuo(im)
c  physical parameters
      parameter(g=g_,rd=rd_,rv=rv_,
     &          cp=cp_,hvap=hvap_,t0c=t0c_)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,
     &          el2orc=hvap*hvap/(rv*cp),eps=rd/rv,epsm1=rd/rv-1.)
c  bounds of parcel origin
      parameter(kliftl=2,kliftu=3)
c  sigma below which to compute moisture convergence
      parameter(sigdhq=0.65)
c  local variables and arrays
c
c  this is checkout demonstration
c
c
c yifeng cui added rsm mp    org:sdsc     date:09-09-05
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
      parameter(ilot=ILOTS,klot=levs_)
      dimension index2(ilot),klcl2(ilot),kbot2(ilot),ktop2(ilot)
      dimension adq2(ilot),aqc2(ilot),atc2(ilot),ps2(ilot),rn2(ilot),
     &          sdq(ilot),sd2(ilot),sdq2(ilot),sdqc2(ilot),
     &          sq2(ilot),sqs2(ilot),sqc2(ilot),stc2(ilot),
     &          dq2(ilot*klot),q2(ilot*klot),qc2(ilot*klot),
     &           t2(ilot*klot),tc2(ilot*klot)
*
#ifdef KUM
      dimension bdq2(ilot), q3(ilot*klot)
#endif
*
c critical relative humidity for computing beta
      rhcrit = 1.
#ifdef KUM
      rhcrit = 1.
#endif
*
c condensation factor to compute sdqc2
      gamma = 1.
#ifdef KUM
      gamma = 1.
#endif
*
c constant beta
#ifdef KUM
      beta0 = 0.
      beta0 = 0.1
#endif
*
c-----------------------------------------------------------------------
c  initialize arrays and compute moisture convergence.
c  compress fields to points with minimum moisture convergence
c  and acceptable boundary layer temperatures.
      do 110 i=1,im
        rn(i)=0.
        kbot(i)=km+1
        ktop(i)=0
        kuo(i)=0
        sdq(i)=0.
110   continue
      do 130 k=1,km
        if(sl(k).ge.sigdhq) then
          do 120 i=1,im
            sdq(i)=sdq(i)+(q(i,k)-qn(i,k))*del(k)
120       continue
        endif
130   continue
c  crconv is in units of g times meters of rain per second.
      crconv=2.e-6
      crconv=crconv*(jcap/80.)**2
      if(jcap.eq.126) crconv=6.e-6
      gmrn=dt*crconv
      im2=0
      do 140 i=1,im
        if(sdq(i)*ps(i).gt.gmrn.and.t(i,2).gt.max(t0c+5.,t(i,3))) then
          im2=im2+1
          index2(im2)=i
        endif
140   continue
      if(im2.eq.0) return
      do 150 i=1,im2
        ps2(i)=ps(index2(i))
150   continue
      do 160 k=1,km
#ifdef CRAY_THREAD
cfpp$ select(vector)
#endif
      do 160 i=1,im2
        ik=(k-1)*im2+i
        dq2(ik)=q(index2(i),k)-qn(index2(i),k)
        q2(ik)=qn(index2(i),k)
        t2(ik)=t(index2(i),k)
160   continue
c-----------------------------------------------------------------------
c  compute moist adiabat and determine cloud boundaries.
c  sum fields within cloud.  restore humidity.
      call mstadb(im2,km,kliftl,kliftu,sl,slk,ps2,t2,q2,
     &            klcl2,kbot2,ktop2,tc2,qc2)
      kbm2=km+1
      ktm2=0
      kbx2=0
      do 210 i=1,im2
        kbm2=min(kbm2,kbot2(i))
        ktm2=max(ktm2,ktop2(i))
c  uncomment next line to turn on evaporation of falling rain
c       if(kbot2(i).le.ktop2(i)) kbx2=max(kbx2,kbot2(i)-1)
        sd2(i)=0.
        sdq2(i)=0.
        sq2(i)=0.
        sqs2(i)=0.
        sqc2(i)=0.
        stc2(i)=0.
#ifdef KUM
        adq2(i) = 0.
        bdq2(i) = 0.
#endif
210   continue
      if(ktm2.lt.kbm2) return
      do 220 k=kbm2,ktm2
      do 220 i=1,im2
        ik=(k-1)*im2+i
        if(k.ge.kbot2(i).and.k.le.ktop2(i)) then
#ifdef ICE
          pvs2 = fpvs(t2(ik))
#else
          pvs2 = fpvs0(t2(ik))
#endif
          qs2=eps*pvs2/(sl(k)*ps2(i)+epsm1*pvs2)
          sd2(i)=sd2(i)+del(k)
          sdq2(i)=sdq2(i)+dq2(ik)*del(k)
          sq2(i)=sq2(i)+q2(ik)*del(k)
          sqs2(i)=sqs2(i)+qs2*del(k)
          sqc2(i)=sqc2(i)+qc2(ik)*del(k)
          stc2(i)=stc2(i)+tc2(ik)*del(k)
        endif
220   continue
      do 230 ik=1,im2*km
#ifdef KUM
        q3(ik) = q2(ik)
#endif
        q2(ik)=q2(ik)+dq2(ik)
230   continue
c-----------------------------------------------------------------------
c  cloud must extend over 0.3 of ps and have moisture convergence.
c  compute partitioning of heating and moistening and rainfall.
c  evaporate rain below cloud base.
c  expand fields back again.
      do 310 i=1,im2
        rn2(i)=0.
        sdqmax=sqc2(i)+cpoel*stc2(i)
        sdqc2(i)=min( gamma*sdq2(i), sdqmax )
        if (sd2(i).gt.0.3.and.sdqc2(i).gt.0..and.sqs2(i).gt.0.) then
          adq2(i)=-sdqc2(i)/sdq2(i)
#ifdef KUM
          if ( adq2(i) .lt. -1. ) then
            adq2(i) = -1.
            bdq2(i) = -(sdqc2(i) - sdq2(i))/sq2(i)
          end if
#endif
          beta = rhcrit - sq2(i)/sqs2(i)
#ifdef KUM
          beta = beta0
#endif
          beta=min(beta,sqc2(i)/sdqc2(i))
          beta=max(beta,1.-cpoel*stc2(i)/sdqc2(i))
          if(beta.le.0.) then
            aqc2(i)=0.
            betanow = 1.
#ifdef KUM
            betanow = 1. - beta
#endif
            atc2(i) = betanow*sdqc2(i)/(cpoel*stc2(i))
#ifdef KUM
            bdq2(i) = bdq2(i) + beta*sdqc2(i)/sq2(i)
#endif
          elseif(beta.ge.1.) then
            aqc2(i)=sdqc2(i)/sqc2(i)
            atc2(i)=0.
          else
            aqc2(i)=beta*sdqc2(i)/sqc2(i)
            atc2(i)=(1.-beta)*sdqc2(i)/(cpoel*stc2(i))
          endif
        else
          kbot2(i)=km+1
          ktop2(i)=0
        endif
310   continue
      do 320 k=kbm2,ktm2
      do 320 i=1,im2
        ik=(k-1)*im2+i
        if(k.ge.kbot2(i).and.k.le.ktop2(i)) then
          dpovg=del(k)/g*ps2(i)
          qchg2=aqc2(i)*qc2(ik)+adq2(i)*dq2(ik)
#ifdef KUM
     &          + bdq2(i)*q3(ik)
#endif
          tchg2=atc2(i)*tc2(ik)
          q2(ik)=q2(ik)+qchg2
          t2(ik)=t2(ik)+tchg2
          rn2(i)=rn2(i)+dpovg*cpoel*tchg2
        endif
320   continue
      do 330 k=kbx2,1,-1
      do 330 i=1,im2
        if(rn2(i).gt.0..and.k.lt.kbot2(i)) then
          ik=(k-1)*im2+i
#ifdef ICE
          pvs2=fpvs(t2(ik))
#else
          pvs2=fpvs0(t2(ik))
#endif
          qs2=eps*pvs2/(sl(k)*ps2(i)+epsm1*pvs2)
          qchg2=qs2-q2(ik)
          if(qchg2.gt.0.) then
            dpovg=del(k)/g*ps2(i)
            qchg2=qchg2/(1.+el2orc*qs2/t2(ik)**2)
            qchg2=qchg2*(1.-exp(-0.32*sqrt(2.*dt*rn2(i))))
            rnchg2=min(dpovg*qchg2,rn2(i))
            qchg2=rnchg2/dpovg
            q2(ik)=q2(ik)+qchg2
            t2(ik)=t2(ik)-elocp*qchg2
            rn2(i)=rn2(i)-rnchg2
          endif
        endif
330   continue
      do 340 i=1,im2
        rn(index2(i))=rn2(i)
        kbot(index2(i))=kbot2(i)
        ktop(index2(i))=ktop2(i)
        if(rn2(i).gt.0.) kuo(index2(i))=1
340   continue
      do 350 k=1,ktm2
#ifdef CRAY_THREAD
cfpp$ select(vector)
#endif
      do 350 i=1,im2
        ik=(k-1)*im2+i
        q(index2(i),k)=q2(ik)
        t(index2(i),k)=t2(ik)
350   continue
c-----------------------------------------------------------------------
      return
      end
