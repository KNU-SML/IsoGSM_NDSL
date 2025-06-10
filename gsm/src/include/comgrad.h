#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#define LATGS latgp_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#define LATGS latg_
#endif
c
      parameter (mcld=3,nseal=2,nbin=100,nlon=2,nlat=4)
      common /comgrad1/ rhcl(nbin,nlon,nlat,mcld,nseal),
     &                  albedr(LONF2S,LATG2S),slmskr(LONF2S,LATG2S),
     &                  coszer(LONF2S,LATG2S),coszdg(LONF2S,LATG2S),
     &                  solc,rsin1,rcos1,rcos2,
     &                  tsmin,tsmax,shmin,shmax,
#ifdef DG
     &                  dtacc,
#endif
     &                  raddt,fjd,r1,alf,dlt,dtlw
      common /comgrad2/ cvr  (LONF2S,LATG2S),cvtr (LONF2S,LATG2S),
     1                  cvbr (LONF2S,LATG2S)
      common /comgrad3/ alvbr(LONF2S,LATG2S),alnbr(LONF2S,LATG2S),
     1                  alvdr(LONF2S,LATG2S),alndr(LONF2S,LATG2S),
#ifndef SW99UPDATE
     2                  paerr(LONF2S,5,LATG2S),
#endif
     3                  avecld(levs_,LATGS),cldl(4,LATGS),
     4                  avecv(3,LATGS),zonht(levs_,LATGS),
     5                  cldsig(levs_,18)
#ifndef SW99UPDATE
      common /comgradi/ istrat,jo3,nalaer,kalb,jsno,itimsw,itimlw,
     1                  icfc,icwp,iswsrc(5),ibnd,ko3
      logical runrad
      common /comgradl/ runrad
      common /alaer1/ alvsf(LONF2S,LATG2S,4),
     1                alnsf(LONF2S,LATG2S,4),
     1                alvwf(LONF2S,LATG2S,4),
     2                alnwf(LONF2S,LATG2S,4),
     2                facsf(LONF2S,LATG2S),
     3                 facwf(LONF2S,LATG2S),
     3                paerf(LONF2S,LATG2S,5)
#else
      common /comgradi/ istrat,jo3,nfalb,nfaer,kalb,jsno,itimsw,itimlw,
     1                  icfc,icwp,iswsrc(5),ibnd,ko3
      logical runrad
      common /comgradl/ runrad
#include <aerparm.h>
      common /alaer1/ alvsf(LONF2S,LATG2S,4),
     1                alnsf(LONF2S,LATG2S,4),
     1                alvwf(LONF2S,LATG2S,4),
     2                alnwf(LONF2S,LATG2S,4),
     2                facsf(LONF2S,LATG2S),  
     3                facwf(LONF2S,LATG2S),
     3                idxcg(nxc,imxae,jmxae),cmixg(nxc,imxae,jmxae),
     4                denng(ndn,imxae,jmxae),kprfg(imxae,jmxae)
#endif
