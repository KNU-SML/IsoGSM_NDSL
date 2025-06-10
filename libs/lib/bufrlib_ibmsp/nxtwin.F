      SUBROUTINE NXTWIN(LUN,IWIN,JWIN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)

      REAL*8 VAL

C----------------------------------------------------------------------
cfpp$ expand (lstrpc)
C----------------------------------------------------------------------

      IF(JWIN.EQ.NVAL(LUN)) THEN
         IWIN = 0
         RETURN
      ENDIF

C  FIND THE NEXT SEQUENTIAL WINDOW
C  -------------------------------

      NODE = INV(IWIN,LUN)
      IF(LSTRPC(NODE,LUN).NE.NODE) print*,'bad node=',node,iwin
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL ABT('NXTWIN - NOT RPC')
      IF(VAL(JWIN,LUN).EQ.0) THEN
         IWIN = 0
      ELSE
         IWIN = JWIN
         JWIN = IWIN+VAL(IWIN,LUN)
      ENDIF

      RETURN
      END
