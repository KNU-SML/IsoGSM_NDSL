#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
c  begin comfiso (for isotopes by kei)
c
      real geshemi1,geshemi2,bengshi1,bengshi2,
     *     dqsfci1,dqsfci2,rresv1,rresv2,
     *     pwati1,pwati2,tmsk,q2mi1,q2mi2
      real geshemc,geshemci1,geshemci2
      real bengshc,bengshci1,bengshci2
      real meqfo,meqfd
      real runoff1,runoff2
      real bgrun1,bgrun2
      common /comfiso/
     *          geshemi1(LONF2S,LATG2S),
     *          geshemi2(LONF2S,LATG2S),
     *          bengshi1(LONF2S,LATG2S),
     *          bengshi2(LONF2S,LATG2S),
     *          dqsfci1(LONF2S,LATG2S),
     *          dqsfci2(LONF2S,LATG2S),
     *          rresv1(LONF2S,LATG2S,lsoil_),
     *          rresv2(LONF2S,LATG2S,lsoil_),
     *          pwati1(LONF2S,LATG2S),
     *          pwati2(LONF2S,LATG2S),
     *          tmsk(LONF2S,LATG2S),
     *          q2mi1(LONF2S,LATG2S),
     *          q2mi2(LONF2S,LATG2S),
     *          geshemc(LONF2S,LATG2S),
     *          geshemci1(LONF2S,LATG2S),
     *          geshemci2(LONF2S,LATG2S),
     *          bengshc(LONF2S,LATG2S),
     *          bengshci1(LONF2S,LATG2S),
     *          bengshci2(LONF2S,LATG2S),
     *          meqfo(LONF2S,LATG2S),
     *          meqfd(LONF2S,LATG2S),
     *          runoff1(LONF2S,LATG2S),
     *          runoff2(LONF2S,LATG2S),
     *          bgrun1(LONF2S,LATG2S),
     *          bgrun2(LONF2S,LATG2S)
c
