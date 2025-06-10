c
c   begin comver
c
      common /rverfidi/
     7   kdt,inistp,limlow,numsum,nummax,ncldb1,
     8   nout0,nout,nrsmi0,nrsmi1,nrsmi2,
     8   nrsmo1,nrsmo2,nrsmop,nrflxi,nrflxf,nrflxp,nrsfcp,
     9   nrflip,nrflop,nrinit,nr2dda,nrpken
      common /rverfid/
     1   am(levs_,levs_),hm(levs_,levs_),tm(levs_,levs_),
     2   bm(levs_,levs_),cm(levs_,levs_),spdmax(levs_),
     3 si(levp1_),sl(levs_),del(levs_),rdel2(levp1_),rmsdot(levm1_),
     4 ci(levp1_),cl(levs_),tov(levs_),   sv(levs_),   rpi(levm1_),
     5 p1(levs_),p2(levs_), h1(levs_),   h2(levs_),rpirec(levm1_),
     6   thour,deltim,dt2,dthr,hdthr,dtltb,dtcvav,dtswav,dtlwav,
     7   rcl,sl1,fhour,shour,dtpost,
     9   filta,filtb,dk,tk,percut,
     x   bndrlx,cowave,dtwave,
chybrid
     x   ak5(levp1_),bk5(levp1_)
#include <rscomcpu.h>
