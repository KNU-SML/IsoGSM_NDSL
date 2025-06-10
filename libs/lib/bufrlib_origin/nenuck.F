      SUBROUTINE NENUCK(NEMO,NUMB,LUN)
#include "ufbprm.h"
      integer TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      CHARACTER*8   NEMO
      CHARACTER*6   NUMB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK TABLE A
C  -------------

      ENTRY NENUAA(NEMO,NUMB,LUN)

      DO N=1,NTBA(LUN)
      IF(NUMB(2:4).EQ.TABA(N,LUN)(1: 3)) GOTO 900
      IF(NEMO     .EQ.TABA(N,LUN)(4:11)) GOTO 900
      ENDDO

      RETURN

C  CHECK TABLE B AND D
C  -------------------

      ENTRY NENUBD(NEMO,NUMB,LUN)

      DO N=1,NTBB(LUN)
      IF(NUMB.EQ.TABB(N,LUN)(1: 6)) GOTO 900
      IF(NEMO.EQ.TABB(N,LUN)(7:14)) GOTO 900
      ENDDO

      DO N=1,NTBD(LUN)
      IF(NUMB.EQ.TABD(N,LUN)(1: 6)) GOTO 900
      IF(NEMO.EQ.TABD(N,LUN)(7:14)) GOTO 900
      ENDDO

      RETURN

C  ERROR EXIT
C  ----------

900   CALL ABORT('NENUCK - DUPLICATE NEM/NUM '//NEMO//' '//NUMB)
      END
