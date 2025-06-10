#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
c  begin comfphys
c
#include <comsfc.h>
c
      real swh,hlw,ozon,solin,sdec,cdec,slag,solhr,clstp
      integer ldum0
      real sfcnsw,sfcdlw
#ifdef VIC
      real sfcusw,sfcdsw
#endif
      real coszen,xlon,xlat,tsflw,dusfc,dvsfc,dtsfc,dqsfc,
     *     dlwsfc,ulwsfc,geshem
#ifdef NOAHYDRO
      real geshem2
#endif
      real dugwd,dvgwd,u10m,v10m,t2m,q2m,psurf,psmean,
     *     bengsh,gflux,slrad,runoff,tmpmax,tmpmin,
     *     ep,cldwrk,hpbl,pwat,snowmelt,snowfall,snowevap,
     *     qull,qvll,alhtfl,evcnp,bgrun,hprime
      common /comfphys/
     *         swh(LONF2S,levs_,LATG2S), hlw(LONF2S,levs_,LATG2S),
     *         ozon(LONF2S,levs_,LATG2S),
     *         solin(LONF2S,levs_,LATG2S), 
     *         sdec,cdec,slag,solhr,clstp,lastep,ldum0,
     *      sfcnsw(LONF2S,LATG2S),    sfcdlw(LONF2S,LATG2S),
#ifdef VIC
     *      sfcusw(LONF2S,LATG2S),    sfcdsw(LONF2S,LATG2S),
#endif
     *      coszen(LONF2S,LATG2S),
     *        xlon(LONF2S,LATG2S),      xlat(LONF2S,LATG2S),
     *        tsflw(LONF2S,LATG2S),
     *       dusfc(LONF2S,LATG2S),     dvsfc(LONF2S,LATG2S),
     *       dtsfc(LONF2S,LATG2S),     dqsfc(LONF2S,LATG2S),
#ifdef CO2PROG
** co2 **
     *     dqsfci1(LONF2S,LATG2S),
#endif
     *      dlwsfc(LONF2S,LATG2S),    ulwsfc(LONF2S,LATG2S),
     *      geshem(LONF2S,LATG2S),
#ifdef NOAHYDRO
     *      geshem2(LONF2S,LATG2S),
#endif
     *       dugwd(LONF2S,LATG2S),     dvgwd(LONF2S,LATG2S),
     *        u10m(LONF2S,LATG2S),      v10m(LONF2S,LATG2S),
     *         t2m(LONF2S,LATG2S),       q2m(LONF2S,LATG2S),
     *       psurf(LONF2S,LATG2S),    psmean(LONF2S,LATG2S),
     *      bengsh(LONF2S,LATG2S),
     *       gflux(LONF2S,LATG2S),     slrad(LONF2S),
     *      runoff(LONF2S,LATG2S),
     *      tmpmax(LONF2S,LATG2S),    tmpmin(LONF2S,LATG2S),
     *          ep(LONF2S,LATG2S),    cldwrk(LONF2S,LATG2S),
     *        hpbl(LONF2S,LATG2S),      pwat(LONF2S,LATG2S),
     *    snowmelt(LONF2S,LATG2S),  snowfall(LONF2S,LATG2S),
     *    snowevap(LONF2S,LATG2S),
     *        qull(LONF2S,LATG2S),      qvll(LONF2S,LATG2S),
     *      alhtfl(LONF2S,LATG2S),     evcnp(LONF2S,LATG2S),
     *      bgrun (LONF2S,LATG2S),
     *      hprime(LONF2S,LATG2S,mtnvar_) 
c
      logical lastep
c
      common /comroms/ 
     *    romsevap(LONF2S,LATG2S), romssens(LONF2S,LATG2S),
     *    romsustr(LONF2S,LATG2S), romsvstr(LONF2S,LATG2S),
     *    romslwup(LONF2S,LATG2S), romsswup(LONF2S,LATG2S),
     *    romslwdn(LONF2S,LATG2S), romsswdn(LONF2S,LATG2S),
     *    romsprcp(LONF2S,LATG2S), romssgz (LONF2S,LATG2S),
     *    romsrtime, romsrtswup, romsrtswdn
