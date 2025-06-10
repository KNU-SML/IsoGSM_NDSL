#define LNWAVP lnwavp
#define LNWAV lnwav
#ifdef MP
#define LNWAVS lnwavp_
#else
#define LNWAVS lnwav_
#endif
c
      common/regbiti/ lab(4), ifin
      common/regbit/ delx, dely, rnnp1max,
     1   snnp1(LNWAVS),rnnp1(LNWAVS),epsx(LNWAVS),epsy(LNWAVS)
      common/regbiti/ icen,igen,icen2,ienst,iensi
      common/regbit/ runid,usrid
      character*8 lab
