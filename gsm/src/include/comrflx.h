c
      real dusfc,dvsfc,dtsfc,dqsfc,dlwsfc,ulwsfc,geshem
#ifdef NOAHYDRO
      real geshem2
#endif
      real u10m,v10m,t2m,q2m,dugwd,dvgwd,psurf,psmean,dtflux
      common /rcomf3/
     *       dusfc(IGRD12S,JGRD12S), dvsfc(IGRD12S,JGRD12S),
     *       dtsfc(IGRD12S,JGRD12S), dqsfc(IGRD12S,JGRD12S),
     *      dlwsfc(IGRD12S,JGRD12S),ulwsfc(IGRD12S,JGRD12S),
     *      geshem(IGRD12S,JGRD12S),
#ifdef NOAHYDRO
     *      geshem2(IGRD12S,JGRD12S),
#endif
     *       u10m(IGRD12S,JGRD12S),v10m(IGRD12S,JGRD12S),
     *       t2m(IGRD12S,JGRD12S),q2m(IGRD12S,JGRD12S),
     *       dugwd(IGRD12S,JGRD12S), dvgwd(IGRD12S,JGRD12S),
     *       psurf(IGRD12S,JGRD12S),psmean(IGRD12S,JGRD12S),
     *       dtflux
c
      real bengsh,gflux,slrad,runoff,tmpmax,tmpmin,
     *     ep,cldwrk,hpbl,pwat,snowmelt,snowfall,snowevap,qull,qvll,
     *     alhtfl,evcnp,bgrun,hprime
      common /rcomf4/
     *        bengsh(IGRD12S,JGRD12S),
     *        gflux(IGRD12S,JGRD12S),  slrad(IGRD12S,JGRD12S),
     *        runoff(IGRD12S,JGRD12S),
     *        tmpmax(IGRD12S,JGRD12S),tmpmin(IGRD12S,JGRD12S),
     *        ep(IGRD12S,JGRD12S),cldwrk(IGRD12S,JGRD12S),
     *        hpbl(IGRD12S,JGRD12S),pwat(IGRD12S,JGRD12S),
     *        snowmelt(IGRD12S,JGRD12S),snowfall(IGRD12S,JGRD12S),
     *        snowevap(IGRD12S,JGRD12S),
     *        qull(IGRD12S,JGRD12S),      qvll(IGRD12S,JGRD12S),
     *        alhtfl(IGRD12S,JGRD12S),evcnp(IGRD12S,JGRD12S),
     *        bgrun (IGRD12S,JGRD12S),
     *        hprime(IGRD12S,JGRD12S,mtnvar_)
c
        common /comroms/
     *        romsevap(LONF2S,LATG2S), romssens(LONF2S,LATG2S),
     *        romsustr(LONF2S,LATG2S), romsvstr(LONF2S,LATG2S),
     *        romslwup(LONF2S,LATG2S), romsswup(LONF2S,LATG2S),
     *        romslwdn(LONF2S,LATG2S), romsswdn(LONF2S,LATG2S),
     *        romsprcp(LONF2S,LATG2S), romssgz (LONF2S,LATG2S),
     *        romsrtime, romsrtswup, romsrtswdn
