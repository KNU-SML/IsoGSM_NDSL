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
      subroutine nfdr_zer
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
#ifdef MP
#include <comrmpi.h>
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#else
#define IGRD12S igrd12_
#define JGRD12S jgrd12_
#endif
#else
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
#include <comnfdr.h>
c
#ifdef ORIGIN_THREAD
c$doacross share(rhmax,rhmin,prcpdur,ppt24,
c$&              t2m2pm,rh2m2pm,cloud2pm,wind2pm,
c$&              tmax24,tmin24)
c$&         local(i,lat)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(rhmax,rhmin,prcpdur,ppt24)
cmic$1 shared(t2m2pm,rh2m2pm,cloud2pm,wind2pm)
cmic$1 shared(tmax24,tmin24)
cmic$1 private(i,lat)
#endif
#ifdef OPENMP
!$omp parallel do private(i,lat)
#endif
      do lat=1,JGRD12S
        do i=1,IGRD12S
          tmax24  (i,lat) = 0.0
          tmin24  (i,lat) = 999.0
          ppt24   (i,lat) = 0.0
          rhmax   (i,lat) = 0.0
          rhmin   (i,lat) = 100.0
          prcpdur (i,lat) = 0.0
          t2m2pm  (i,lat) = -999.
          rh2m2pm (i,lat) = -999.
          cloud2pm(i,lat) = -999.
          wind2pm (i,lat) = -999.
        enddo
      enddo
c
      return
      end
