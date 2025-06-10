#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
      common /radiag/ fluxr(LONF2S,LATG2S,26)
      common /radiag/ cvavg(LONF2S,LATG2S)
