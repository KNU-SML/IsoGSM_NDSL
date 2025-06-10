      SUBROUTINE UFBCUP(LUBIN,LUBOT)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
 
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
 
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
 
      CHARACTER*10 TAG,TAGI(INV_),TAGO
      CHARACTER*3  TYP
      DIMENSION    NINI(INV_)
      REAL*8       VAL
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C  CHECK THE FILE STATUSES AND I-NODE
C  ----------------------------------
 
      CALL STATUS(LUBIN,LUI,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IM.EQ.0) GOTO 901
      IF(INODE(LUI).NE.INV(1,LUI)) GOTO 902
 
      CALL STATUS(LUBOT,LUO,IL,IM)
      IF(IL.EQ.0) GOTO 903
      IF(IM.EQ.0) GOTO 904
 
C  MAKE A LIST OF UNIQUE TAGS IN INPUT BUFFER
C  ------------------------------------------
 
      NTAG = 0
 
      DO 5 NI=1,NVAL(LUI)
      NIN = INV(NI,LUI)
      IF(ITP(NIN).GE.2) THEN
         DO NV=1,NTAG
         IF(TAGI(NV).EQ.TAG(NIN)) GOTO 5
         ENDDO
         NTAG = NTAG+1
         NINI(NTAG) = NI
         TAGI(NTAG) = TAG(NIN)
      ENDIF
5     ENDDO
 
      IF(NTAG.EQ.0) GOTO 905
 
C  GIVEN A list MAKE ONE COPY OF COMMON ELEMENTS TO OUTPUT BUFFER
C  --------------------------------------------------------------
 
      DO 10 NV=1,NTAG
      NI = NINI(NV)
      DO NO=1,NVAL(LUO)
      TAGO = TAG(INV(NO,LUO))
      IF(TAGI(NV).EQ.TAGO) THEN
         VAL(NO,LUO) = VAL(NI,LUI)
         GOTO 10
      ENDIF
      ENDDO
10    ENDDO
 
C  ALL EXITS HERE
C  --------------
 
      RETURN
900   CALL ABORT('UFBCUP - INPUT  FILE IS NOT OPEN             ')
901   CALL ABORT('UFBCUP - INPUT  MESG IS NOT OPEN             ')
902   CALL ABORT('UFBCUP - INPUT  I-NODE  MISMATCH             ')
903   CALL ABORT('UFBCUP - OUTPUT FILE IS NOT OPEN             ')
904   CALL ABORT('UFBCUP - OUTPUT MESG IS NOT OPEN             ')
905   CALL ABORT('UFBCUP - NO TAGS IN INPUT BUFFER             ')
      END
