#ifdef MP
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#else
#define IGRD12S igrd12_ 
#define JGRD12S jgrd12_
#endif 
c
      parameter (mcld=3,nseal=2,nbin=100,nlon=2,nlat=4)
      common /comgrad/ rhcl(nbin,nlon,nlat,mcld,nseal),
     &                 solc,rsin1,rcos1,rcos2,
     &                 tsmin,tsmax,shmin,shmax,
     &                 raddt,fjd,r1,alf,dlt,dlon,dtlw,rrs2
      common /comgrad/ cvr (IGRD12S,JGRD12S),
     1                  cvtr(IGRD12S,JGRD12S),
     1                 cvbr(IGRD12S,JGRD12S)
      common /comgrad/ alvbr(IGRD12S,JGRD12S),
     1                 alnbr(IGRD12S,JGRD12S),
     1                 alvdr(IGRD12S,JGRD12S),
     2                 alndr(IGRD12S,JGRD12S),
#ifdef RSMVECTORIZE
     2                 paerr(igrd12p_*jgrd12p_,5,1),
#else
     2                 paerr(IGRD12S,5,JGRD12S),
#endif
     3                 avecld(levs_,jgrd1_),cldl(4,jgrd1_),
     4                 avecv(3,jgrd1_),zonht(levs_,jgrd1_),
     5                 cldsig(levs_,18)
      common /alaer1/alvsf(IGRD12S,JGRD12S,4),
     1               alnsf(IGRD12S,JGRD12S,4),
     1               alvwf(IGRD12S,JGRD12S,4),
     2               alnwf(IGRD12S,JGRD12S,4),
     2               facsf(IGRD12S,JGRD12S),
     3               facwf(IGRD12S,JGRD12S),
     3               paerf(IGRD12S,JGRD12S,5)
      common /comgradi/idtln,idtls,istrat,jo3,nalaer,
     &                 kalb,jsno,itimsw,itimlw,
     1                 icfc,icwp,iswsrc(5),ibnd,ko3
      logical runrad
      common /comgradl/runrad
