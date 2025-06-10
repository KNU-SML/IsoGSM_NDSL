      FUNCTION LSTRPS(NODE,LUN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*10 TAG
      CHARACTER*3  TYP

C----------------------------------------------------------------------
C----------------------------------------------------------------------

      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) GOTO 900

      NOD = NODE

C  FIND THIS OR THE PREVIOUS RPS NODE
C  ----------------------------------

10    IF(TYP(NOD).NE.'RPS') THEN
         NOD = JMPB(NOD)
         IF(NOD.NE.0) GOTO 10
      ENDIF

      LSTRPS = NOD

      RETURN
900   CALL ABORT('LSTRPS - TABLE NODE IS OUT OF BOUNDS')
      END
