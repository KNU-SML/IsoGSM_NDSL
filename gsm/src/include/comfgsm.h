c
      common /comgsmi/ lnts2,lnoffset,lons2,lats2,n1,
#ifdef RSM
     3 nrsmi1,nrsmi2,nrflip,
     & nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken,
#endif
#ifdef LFM
     & nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,
#endif
     5 kpfix,ksfcx,ksig,ksfc,kpost,krestart,klfm,krsm
#ifdef COUPLE
     6,kcouple
#endif
      common /comgsmr/ solsec,
#ifdef LFM
     & weight,filtwin,
#endif
     2 dthr,hdthr
