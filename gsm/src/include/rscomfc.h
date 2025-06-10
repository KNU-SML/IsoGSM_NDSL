#define LNUVS lnuv_
#define LNT22S lnt22_
#define CLNGRDS clngrd_ 
#include <machine.h>
c
c  begin (rscomfc)
c    version with stacked transforms
c
#ifdef G2R
       common /rcomfgi/ idate(4)
       common /rcomfg/ si(levs_+1),sl(levs_),
     *             eps(LNUVS),epsi(LNUVS),
     *              gz(LNT22S),
     *              ze(LNT22S,levs_),
     *              di(LNT22S,levs_),
     *              te(LNT22S,levs_),
     *              rq(LNT22S,levh_),
     *               q(LNT22S)
#endif
#ifdef C2R
       common /rcomfci/ idate(4)
       common /rcomfc/ si(levs_+1),sl(levs_)
     1 ,  q(CLNGRDS)
     2 , te(CLNGRDS,levs_)
     3 , rq(CLNGRDS,levh_)
     4 , gz(CLNGRDS)
     5 , uu(CLNGRDS,levs_)
     6 , vv(CLNGRDS,levs_)
     7,flat(CLNGRDS), flon(CLNGRDS)
     8, fm2(CLNGRDS), fm2x(CLNGRDS), fm2y(CLNGRDS)
#endif
c
      common /rcomfc/ filta
#ifdef CRA
      common /rcomfci/ ncpus,ncpus1
#endif
