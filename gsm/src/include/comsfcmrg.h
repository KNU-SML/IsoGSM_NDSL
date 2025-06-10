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
c  this include file is to access sfc variables as one array
c
      real sfcfcs
      common /comsfc/sfcfcs(LONF2S*LATG2S,numsfcs)
