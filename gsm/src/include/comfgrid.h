#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
c  begin comfgrid
c
       common /comfgrid/
     * colrad(latg2_),wgt(latg2_),wgtcs(latg2_),rcs2(latg2_),
     * sinlat(latg2_),
     * colrab(LATG2S),wgb(LATG2S),wgbcs(LATG2S),rbs2(LATG2S),
     * sinlab(LONF2S,LATG2S),coslab(LONF2S,LATG2S)
#ifndef MP
     *,lwvdef(jcap_+1),latdef(latg_/2)
#endif
