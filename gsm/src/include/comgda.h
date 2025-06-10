#ifdef DG3
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
c  diagnostic indexes and flags
c
      parameter(kdgda=13)
      parameter(kdtlarg=1,kdtconv=2,kdqconv=3,kdtshal=4,kdqshal=5,
     &          kdtvrdf=6,kduvrdf=7,kdvvrdf=8,kdqvrdf=9,
     &          kdthsw=10,kdthlw=11,kdtcld=12,kdtccv=13)
      character*8 cnmgda
      common /comgdc/ cnmgda(kdgda)
      parameter(ntgda=92)
#ifdef RSMVECTORIZE
      parameter(nwgda=((igrd12p_*jgrd12p_*levs_-1)/512+1)*512)
#else
      parameter(nwgda=((LONF2S*levs_-1)/512+1)*512)
#endif
      parameter(nrgda=LATG2S)
      common /comgda/ ipugda(kdgda),ibmgda(kdgda)
      common /comgdd/ gdd(nwgda*kdgda*nrgda)
#endif
