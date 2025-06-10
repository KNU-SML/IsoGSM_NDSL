      SUBROUTINE NEMTBA(LUN,NEMO,MTYP,MSBT,INOD)
#include "ufbprm.h"
      integer TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)

      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      CHARACTER*20  NEMT

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      NEMT = NEMO
      IRET = 0

C  LOOK FOR NEMO IN TABLE A
C  ------------------------

      DO I=1,NTBA(LUN)
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN
         MTYP = IDNA(I,LUN,1)
         MSBT = IDNA(I,LUN,2)
         INOD = MTAB(I,LUN)
         IF(MTYP.LT.0 .OR. MTYP.GT.255) GOTO 900
         IF(MSBT.LT.0 .OR. MSBT.GT.255) GOTO 901
         RETURN
      ENDIF
      ENDDO

      CALL ABORT('NEMTBA - CANT FIND '//NEMT)
900   CALL ABORT('NEMTBA - BAD MTYP  '//NEMT)
901   CALL ABORT('NEMTBA - BAD MSBT  '//NEMT)
      END
