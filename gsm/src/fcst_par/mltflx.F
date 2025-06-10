#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine mltflx(fac,dusfc,dvsfc,dtsfc,dqsfc,dlwsfc,ulwsfc,
     1 bengsh,gflux,
     2 dugwd,dvgwd,psmean,
     3 dtfulx)
#include <paramodel.h>
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
      dimension dusfc(LONF2S,LATG2S),dvsfc(LONF2S,LATG2S)
      dimension dtsfc(LONF2S,LATG2S),dqsfc(LONF2S,LATG2S)
      dimension dlwsfc(LONF2S,LATG2S),ulwsfc(LONF2S,LATG2S)
      dimension bengsh(LONF2S,LATG2S),gflux(LONF2S,LATG2S)
      dimension dugwd(LONF2S,LATG2S),dvgwd(LONF2S,LATG2S)
      dimension psmean(LONF2S,LATG2S)
      do 23500 l=1,LATG2S
      do 23500 j=1,LONF2S
       dusfc(j,l)=fac*dusfc(j,l)
       dvsfc(j,l)=fac*dvsfc(j,l)
       dtsfc(j,l)=fac*dtsfc(j,l)
       dqsfc(j,l)=fac*dqsfc(j,l)
      dlwsfc(j,l)=fac*dlwsfc(j,l)
      ulwsfc(j,l)=fac*ulwsfc(j,l)
      bengsh(j,l)=fac*bengsh(j,l)
       gflux(j,l)=fac*gflux(j,l)
       dugwd(j,l)=fac*dugwd(j,l)
       dvgwd(j,l)=fac*dvgwd(j,l)
       psmean(j,l)=fac*psmean(j,l)
23500 continue
c--------------------------------------------------------------
      return
      end
