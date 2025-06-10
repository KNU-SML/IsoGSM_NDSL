      FUNCTION LSTJPB(NODE,LUN,JBTYP)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*(*) JBTYP
      CHARACTER*10 TAG
      CHARACTER*3  TYP

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C  MAKE SURE WE ARE ALL ON THE SAME PAGE
C  -------------------------------------

      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) THEN
         PRINT*,INODE(LUN),':',NODE,':',TAG(NODE)
         CALL ABT('LSTJPB - TABLE NODE IS OUT OF BOUNDS')
      ENDIF

C  FIND THIS OR THE PREVIOUS RPC NODE
C  ----------------------------------

      LSTJPB = NODE

10    IF(TYP(LSTJPB).NE.JBTYP) THEN
         LSTJPB = JMPB(LSTJPB)
         IF(LSTJPB.NE.0) GOTO 10
      ENDIF

      RETURN
      END
