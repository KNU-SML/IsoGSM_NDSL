c
c  begin comriso (for isotopes by kei)
c
      real geshemi1,geshemi2,bengshi1,bengshi2,
     *     dqsfci1,dqsfci2,rresv1,rresv2,
     *     pwati1,pwati2,tmsk,q2mi1,q2mi2
      real geshemc,geshemci1,geshemci2
      real bengshc,bengshci1,bengshci2
      real meqfo,meqfd
      real runoff1,runoff2
      real bgrun1,bgrun2
      common /comriso/
     *          geshemi1(IGRD12S,JGRD12S),
     *          geshemi2(IGRD12S,JGRD12S),
     *          bengshi1(IGRD12S,JGRD12S),
     *          bengshi2(IGRD12S,JGRD12S),
     *          dqsfci1(IGRD12S,JGRD12S),
     *          dqsfci2(IGRD12S,JGRD12S),
     *          rresv1(IGRD12S,JGRD12S,lsoil_),
     *          rresv2(IGRD12S,JGRD12S,lsoil_),
     *          pwati1(IGRD12S,JGRD12S),
     *          pwati2(IGRD12S,JGRD12S),
     *          tmsk(IGRD12S,JGRD12S),
     *          q2mi1(IGRD12S,JGRD12S),
     *          q2mi2(IGRD12S,JGRD12S),
     *          geshemc(IGRD12S,JGRD12S),
     *          geshemci1(IGRD12S,JGRD12S),
     *          geshemci2(IGRD12S,JGRD12S),
     *          bengshc(IGRD12S,JGRD12S),
     *          bengshci1(IGRD12S,JGRD12S),
     *          bengshci2(IGRD12S,JGRD12S),
     *          meqfo(IGRD12S,JGRD12S),
     *          meqfd(IGRD12S,JGRD12S),
     *          runoff1(IGRD12S,JGRD12S),
     *          runoff2(IGRD12S,JGRD12S),
     *          bgrun1(IGRD12S,JGRD12S),
     *          bgrun2(IGRD12S,JGRD12S)
c
