#include <machine.h>
c
c  begin comfver
c
      common/comfver/am(levs_,levs_),hm(levs_,levs_),tm(levs_,levs_),
     1               bm(levs_,levs_),cm(levs_,levs_),spdmax(levs_),
     2 si(levp1_),sl(levs_),del(levs_),rdel2(levs_),rmsdot(levm1_),
     3 ci(levp1_),cl(levs_),tov(levs_),gv(levs_),sv(levs_),rpi(levm1_),
     4 p1(levs_),p2(levs_), h1(levs_),   h2(levs_),rpirec(levm1_),
     5   thour,deltim,sl1,z00,fhour,shour,dtcvav,dtpost,
     6   filta,filtb,dk,tk,percut,dtswav,dtlwav,avprs0,
     7   cowave,dtwave
chybrid
     8  ,ak5(levp1_),bk5(levp1_)
      common/comfveri/ kdt,jdt,inistp,limlow,maxstp,numsum,nummax
#ifdef DYNAMIC_ALLOC
     8   ,ncldb1,ncpus,ncpus1
#endif
      common/comfverl/ stepone
      logical stepone
