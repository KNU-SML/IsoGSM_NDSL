#include <rsparltb.h>
#ifdef MP
#define LNGRDS lngrdp_
#else
#define LNGRDS lngrd_
#endif
c
c  begin comltb
c
      common/rltbcom/ 
     1   rltb(LNGRDS),  sltb(LNGRDS)
      common/ltbcom/ 
     1   blat(lngrdb), blon(lngrdb)
     1  ,bcsln(lngrdb),bsnln(lngrdb)
     1  , wsltb(bgf_,1-border_:border_)
     1  ,dwsltb(bgf_,1-border_:border_)
     2  ,             qb(lngrdb)
     2  ,uub(lngrdb,levs_),vvb(lngrdb,levs_)
     2  ,teb(lngrdb,levs_),rqb(lngrdb,levh_)
     2  , gzb(lngrdb)
     3  ,             qltb(lngrdb)
     3  ,uultb(lngrdb,levs_),vvltb(lngrdb,levs_)
     3  ,teltb(lngrdb,levs_),rqltb(lngrdb,levh_)
