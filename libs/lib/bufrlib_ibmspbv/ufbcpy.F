      SUBROUTINE UFBCPY(LUBIN,LUBOT)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)

      REAL*8 VAL

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
      IF(INODE(LUI).NE.INODE(LUO)) GOTO 905

C  EVERYTHING OKAY COPY USER ARRAY FROM LUI TO LUO
C  -----------------------------------------------

      NVAL(LUO) = NVAL(LUI)

      DO N=1,NVAL(LUI)
      INV(N,LUO) = INV(N,LUI)
      VAL(N,LUO) = VAL(N,LUI)
      ENDDO

      RETURN
900   CALL ABT('UFBCPY - INPUT  FILE IS NOT OPEN             ')
901   CALL ABT('UFBCPY - INPUT  MESG IS NOT OPEN             ')
902   CALL ABT('UFBCPY - INPUT  I-NODE  MISMATCH             ')
903   CALL ABT('UFBCPY - OUTPUT FILE IS NOT OPEN             ')
904   CALL ABT('UFBCPY - OUTPUT MESG IS NOT OPEN             ')
905   CALL ABT('UFBCPY - IN/OUT I-NODE  MISMATCH             ')
      END
