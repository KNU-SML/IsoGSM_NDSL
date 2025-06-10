      SUBROUTINE GETWIN(NODE,LUN,IWIN,JWIN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)

      REAL*8 VAL

C----------------------------------------------------------------------
cfpp$ expand (lstrpc)
C----------------------------------------------------------------------

      IRPC = LSTRPC(NODE,LUN)

      IF(IRPC.EQ.0) THEN
         IWIN = INVWIN(NODE,LUN,JWIN,NVAL(LUN))
         IF(IWIN.EQ.0 .and. jwin.gt.1) RETURN
         IWIN = 1
         JWIN = NVAL(LUN)
         RETURN
      ELSE
         IWIN = INVWIN(IRPC,LUN,JWIN,NVAL(LUN))
         IF(IWIN.EQ.0) THEN
            RETURN
         ELSEIF(VAL(IWIN,LUN).EQ.0.) THEN
            IWIN = 0
            RETURN
         ENDIF
      ENDIF

      JWIN = INVWIN(IRPC,LUN,IWIN+1,NVAL(LUN))
      IF(JWIN.EQ.0) CALL ABT('GETWIN - MISSING BRACKET')

      RETURN
      END
