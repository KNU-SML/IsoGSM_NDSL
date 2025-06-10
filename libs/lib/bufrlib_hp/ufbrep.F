      SUBROUTINE UFBREP(LUNIO,USR,I1,I2,IRET,STR)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
 
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
 
      CHARACTER*(*) STR
      DIMENSION     USR(I1,I2)
      REAL*8        USR,VAL
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C  CHECK THE FILE STATUS AND I-NODE
C  --------------------------------
 
      LUNIT = ABS(LUNIO)
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IM.EQ.0) GOTO 901
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902
      IO = MIN(MAX(0,IL),1)
      IF(LUNIO.NE.LUNIT) IO = 0
 
C  PARSE THE INPUT STRING - READ/WRITE VALUES
C  ------------------------------------------
 
      CALL STRING(STR,LUN,I1,IO)
      CALL UFBRP(LUN,USR,I1,I2,IO,IRET)
 
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN
         PRINT*,STR
         GOTO 903
      ENDIF
 
      RETURN
900   CALL ABORT('UFBREP - FILE IS CLOSED                     ')
901   CALL ABORT('UFBREP - NO MESSAGE OPEN                    ')
902   CALL ABORT('UFBREP - I-NODE MISMATCH                    ')
903   CALL ABORT('UFBREP - INCOMPLETE WRITE                   ')
      END
