#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
c
c  begin comfspec
c
       common /comfspec/ idate(4),relvor(levs_),absvor(levs_),
     1 eps(lnuv_),epsi(lnuv_),gz(LNT22S)
       common /comfspec/
     *    zem    (LNT22S,levs_),
     *    dim    (LNT22S,levs_),
     *    tem    (LNT22S,levs_),
#ifndef NISLQ
     *    rm     (LNT22S,levh_),
#endif
     *    qm     (LNT22S)
       common /comfspec/
     *    ze     (LNT22S,levs_),
     *    di     (LNT22S,levs_),
     *    te     (LNT22S,levs_),
#ifndef NISLQ
     *    rq     (LNT22S,levh_),
#endif
     *    uln    (LNT22S,levs_),
     *    vln    (LNT22S,levs_),
     *    dpdphi (LNT22S),
     *    dpdlam (LNT22S),
     *    q      (LNT22S),
     *    qlap   (LNT22S)
       common /comfspec/
     *    z      (LNT22S),
     *    uu     (LNT22S,levs_),
     *    vv     (LNT22S,levs_),
     *    y      (LNT22S,levs_),
#ifndef NISLQ
     *    rt     (LNT22S,levh_),
#endif     
     *    x      (LNT22S,levs_),
     *    w      (LNT22S,levs_)
#ifdef NISLQ
       common /comfspec/ rm (LNT22S,levh_)
       common /comfspec/ rq (LNT22S,levh_)
       common /comfspec/ rt (LNT22S,levh_)
#endif
