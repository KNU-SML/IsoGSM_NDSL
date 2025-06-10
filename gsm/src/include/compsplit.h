c
#ifndef RSMVECTORIZE
      common/compsplit/ dudtm(LONF22S,levs_,LATG2S),
     *                  dvdtm(LONF22S,levs_,LATG2S),
     *                  dtdtm(LONF22S,levs_,LATG2S),
     *                  drdtm(LONF22S,levh_,LATG2S),
     *                  dpdtm(LONF22S,      LATG2S)
#else
      common/compsplit/ dudtm(igrd12p_*jgrd12p_,levs_,1),
     *                  dvdtm(igrd12p_*jgrd12p_,levs_,1),
     *                  dtdtm(igrd12p_*jgrd12p_,levs_,1),
     *                  drdtm(igrd12p_*jgrd12p_,levh_,1),
     *                  dpdtm(igrd12p_*jgrd12p_,      1)
#endif
c
