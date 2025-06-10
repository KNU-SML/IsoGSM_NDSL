      FUNCTION INVTAG(NODE,LUN,INV1,INV2)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(TAB_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
 
      CHARACTER*10 TAG,tagn
      CHARACTER*3  TYP
      REAL*8       VAL
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      INVtag = 0
      IF(NODE.EQ.0) RETURN
      tagn = tag(node)
 
C  SEARCH BETWEEN INV1 AND INV2
C  ----------------------------
 
10    DO INVtag=INV1,INV2
      IF(tag(INV(INVtag,LUN)).EQ.tagn) RETURN
      ENDDO
 
      INVtag = 0
      RETURN
      END
