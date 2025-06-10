      SUBROUTINE CONWIN(LUN,INC1,INC2,NBMP)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10)

      REAL*8 VAL
      REAL*4 VALS

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C  SPECIAL CASES
C  -------------

      IF(NCON.EQ.0) THEN
         INC1 = 1
         INC2 = NVAL(LUN)
         RETURN
      ENDIF

      IF(INC1.GT.1 .AND. KONS(NCON).EQ.5) THEN
         CALL NXTWIN(LUN,INC1,INC2)
         RETURN
      ENDIF

C  EVALUATE CONDITIONS TO SEE IF ANY MORE CASES
C  --------------------------------------------

10    DO NC=1,NCON
      IF(KONS(NC).EQ.5) THEN
         INC1 = INVWIN(NODC(NC),LUN,INC1,NVAL(LUN))
         CALL USRTPL(LUN,INC1-1,NBMP)
         CALL NEWWIN(LUN,INC1,INC2)
      ELSE
15       CALL GETWIN(NODC(NC),LUN,INC1,INC2)
         IF(INC1.EQ.0 .AND. NC.EQ.1) RETURN
         IF(INC1.EQ.0              ) GOTO10
         ICON = INVCON(NC,LUN,INC1,INC2)
         IF(ICON.EQ.0) GOTO 15
      ENDIF
      ENDDO

      RETURN
      END
