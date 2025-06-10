#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#include <define.h>
      subroutine zerflx(ifluxr)
c
#include <paramodel.h>
c
c regional spectral model developed by hann-ming henry juang
c       version with nested to global spectral model
c
c. this routine setup all the routine and constant for rsm.
c    and get input data from gsm or others include surface files.
c
#ifdef RSM
#include <rscomf.rerun.h>
#ifdef MP
#include <comrmpi.h>
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#else
#define IGRD12S igrd12_
#define JGRD12S jgrd12_
#endif
#else
#include <comfphys.h>
#ifdef ISOTOPE
#include <comfiso.h>
#endif
#include <radiag.h>
#include <comfver.h>
#ifdef MP
#include <commpi.h>
#define IGRD12S lonf2p_
#define JGRD12S latg2p_
#else
#define IGRD12S lonf2_
#define JGRD12S latg2_
#endif
#endif
c
#ifndef RSM
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#else
#ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(dusfc,dvsfc,dtsfc,dqsfc,
c$&              geshem,bengsh,ulwsfc,dlwsfc,
c$&              gflux,tmpmax,tmpmin,runoff,ep,
c$&              cldwrk,dugwd,dvgwd,psmean,cvavg,
c$&              snowmelt,snowevap,snowfall,qull,qvll,
c$&              alhtfl,evcnp,bgrun
c$&         local(i,lat)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(dusfc,dvsfc,dtsfc,dqsfc)
cmic$1 shared(geshem,bengsh,ulwsfc,dlwsfc)
cmic$1 shared(gflux,tmpmax,tmpmin,runoff,ep)
cmic$1 shared(cldwrk,dugwd,dvgwd,psmean,cvavg)
cmic$1 shared(snowmelt,snowevap,snowfall)
cmic$1 shared(qull,qvll,alhtfl,evcnp,bgrun)
cmic$1 private(i,lat)
#endif
#ifdef OPENMP
!$omp parallel do private(i,lat)
#endif
c
      dhour=dtpost/3600.
      if(dtpost.gt.0) then
        rtime=1./dtpost
      else
        rtime=0.
      endif
c
      secswr=max(dhour,dtswav) * 3600.
      seclwr=max(dhour,dtlwav) * 3600.
      if(secswr.gt.0.) then
        rtimsw=1./secswr
      else
        rtimsw=1.
      endif
      if(seclwr.gt.0.) then
        rtimlw=1./seclwr
      else
        rtimlw=1.
      endif
c
      if(ifluxr.eq.3) then
       do lat=1,JGRD12S
        do i=1,IGRD12S
         romsevap(i,lat)=romsevap(i,lat)+dqsfc(i,lat)
         romssens(i,lat)=romssens(i,lat)+dtsfc(i,lat)
         romsustr(i,lat)=romsustr(i,lat)+dusfc(i,lat)
         romsvstr(i,lat)=romsvstr(i,lat)+dvsfc(i,lat)
         romslwup(i,lat)=romslwup(i,lat)+ulwsfc(i,lat)
         romsswup(i,lat)=romsswup(i,lat)+fluxr(i,lat,3)
         romslwdn(i,lat)=romslwdn(i,lat)+dlwsfc(i,lat)
         romsswdn(i,lat)=romsswdn(i,lat)+fluxr(i,lat,4)
         romsprcp(i,lat)=romsprcp(i,lat)+geshem(i,lat)
        enddo
       enddo
       romsrtime =romsrtime + 1./rtime
       romsrtswup=romsrtswup+ 1./rtimsw
       romsrtswdn=romsrtswdn+ 1./rtimsw
      endif
c
      if(ifluxr.ne.2) then
       do lat=1,JGRD12S
        do i=1,IGRD12S
          geshem(i,lat) = 0.0 e 0
          bengsh(i,lat) = 0.0 e 0
          dusfc (i,lat) = 0.0 e 0
          dvsfc (i,lat) = 0.0 e 0
          dtsfc (i,lat) = 0.0 e 0
          dqsfc (i,lat) = 0.0 e 0
          ulwsfc(i,lat) = 0.0 e 0
          dlwsfc(i,lat) = 0.0 e 0
          gflux (i,lat) = 0.0 e 0
          tmpmax(i,lat) = 0.0 e 0
          tmpmin(i,lat) = 1.0 e 10
          runoff(i,lat) = 0.0 e 0
          ep    (i,lat) = 0.0 e 0
          cldwrk(i,lat) = 0.0 e 0
          dugwd (i,lat) = 0.0 e 0
          dvgwd (i,lat) = 0.0 e 0
          cvavg (i,lat) = 0.0 e 0
          psmean(i,lat) = 0.0 e 0
          snowmelt(i,lat) = 0.0 e 0
          snowevap(i,lat) = 0.0 e 0
          snowfall(i,lat) = 0.0 e 0
          qull(i,lat) = 0.0 e 0
          qvll(i,lat) = 0.0 e 0
          alhtfl(i,lat) = 0.0 e 0
          evcnp(i,lat) = 0.0 e 0
          bgrun(i,lat) = 0.0 e 0
#ifdef NOAHYDRO
          geshem2(i,lat) = 0.0 e 0
#endif
#ifdef CO2PROG
** co2 **
          dqsfci1(i,lat) = 0.0 e 0
#endif
#ifdef ISOTOPE
*** isotope ***
          geshemi1(i,lat) = 0.0 e 0
          geshemi2(i,lat) = 0.0 e 0
          bengshi1(i,lat) = 0.0 e 0
          bengshi2(i,lat) = 0.0 e 0
          dqsfci1(i,lat) = 0.0 e 0
          dqsfci2(i,lat) = 0.0 e 0
          geshemc(i,lat) = 0.0 e 0
          geshemci1(i,lat) = 0.0 e 0
          geshemci2(i,lat) = 0.0 e 0
          bengshc(i,lat) = 0.0 e 0
          bengshci1(i,lat) = 0.0 e 0
          bengshci2(i,lat) = 0.0 e 0
          meqfo(i,lat)= 1.0 
          meqfd(i,lat)= 1.0 
          runoff1(i,lat) = 0.0 e 0
          runoff2(i,lat) = 0.0 e 0
          bgrun1(i,lat) = 0.0 e 0
          bgrun2(i,lat) = 0.0 e 0
#endif
        enddo
       enddo
       dtflux=0.0 e 0
      endif
c
      if((ifluxr.eq.0).or.(ifluxr.eq.3)) then
      do iv=1,26
        do lat=1,JGRD12S
          do i=1,IGRD12S
            fluxr(i,lat,iv) = 0.e0
          enddo
        enddo
      enddo
      do lat=1,JGRD12S
        do i=1,IGRD12S
          dlwsfc(i,lat)=0.
          ulwsfc(i,lat)=0.
          cvavg (i,lat)=0.
        enddo
      enddo
      endif
c
      if(ifluxr.ne.3) then
       do lat=1,JGRD12S
        do i=1,IGRD12S
          romsevap(i,lat)=0.
          romssens(i,lat)=0.
          romsustr(i,lat)=0.
          romsvstr(i,lat)=0.
          romslwup(i,lat)=0.
          romsswup(i,lat)=0.
          romslwdn(i,lat)=0.
          romsswdn(i,lat)=0.
          romsprcp(i,lat)=0.
        enddo
       enddo
       romsrtime=0.
       romsrtswup=0.
       romsrtswdn=0.
      endif
c
      return
      end
