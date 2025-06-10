      SUBROUTINE NEWWIN(LUN,IWIN,JWIN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)

      REAL*8 VAL

C----------------------------------------------------------------------
cfpp$ expand (lstrpc)
C----------------------------------------------------------------------

      IF(IWIN.EQ.1) THEN
         JWIN = NVAL(LUN)
         RETURN
      ENDIF

C  REFIND THE JWIN BOUNDARY FROM IWIN
C  ----------------------------------

      NODE = INV(IWIN,LUN)
      IF(LSTRPC(NODE,LUN).NE.NODE) CALL ABT('NEWWIN - NOT RPC')
      JWIN = IWIN+VAL(IWIN,LUN)

      RETURN
      END
