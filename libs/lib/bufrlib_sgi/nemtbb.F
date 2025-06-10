      SUBROUTINE NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)
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
      CHARACTER*24  UNIT
      CHARACTER*8   NEMO

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF(ITAB.LE.0 .OR. ITAB.GT.NTBB(LUN)) GOTO 900

C  PULL OUT TABLE B INFORMATION
C  ----------------------------

      IDN  = IDNB(ITAB,LUN)
      NEMO = TABB(ITAB,LUN)( 7:14)
      UNIT = TABB(ITAB,LUN)(71:94)
      ISCL = VALS(TABB(ITAB,LUN)( 95: 98))
      IREF = VALS(TABB(ITAB,LUN)( 99:109))
      IBIT = VALS(TABB(ITAB,LUN)(110:112))

C  CHECK TABLE B CONTENTS
C  ----------------------

      IF(IDN.LT.IFXY('000000')) GOTO 901
      IF(IDN.GT.IFXY('063255')) GOTO 901

      IF(ISCL.LT.-999 .OR. ISCL.GT.999) GOTO 902
      IF(IREF.LE.-1E9 .OR. IREF.GE.1E9) GOTO 903
      IF(IBIT.LE.   0 .OR. IBIT.GT. 64) GOTO 904

      IF(UNIT.EQ.'CCITT IA5' .AND. MOD(IBIT,8).NE.0) GOTO 905

      RETURN
900   CALL ABORT('NEMTBB - ITAB NOT IN TABLE B'         )
901   CALL ABORT('NEMTBB - BAD DESCRIPTOR VALUE: '//NEMO)
902   CALL ABORT('NEMTBB - BAD SCALE VALUE     : '//NEMO)
903   CALL ABORT('NEMTBB - BAD REFERENCE VALUE : '//NEMO)
904   CALL ABORT('NEMTBB - BAD BIT WIDTH       : '//NEMO)
905   CALL ABORT('NEMTBB - BAD CHAR BIT WIDTH  : '//NEMO)
      END
