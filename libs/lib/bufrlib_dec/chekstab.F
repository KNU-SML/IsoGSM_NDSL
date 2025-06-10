      SUBROUTINE CHEKSTAB(LUN)
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
      CHARACTER*8   NEMO,NEMS(TBD_)
      CHARACTER*1   TAB
      DIMENSION     IRPS(TBD_),KNTS(TBD_)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  THERE MUST BE ENTRIES IN TABLES A, B, AND D
C  -------------------------------------------

      IF(NTBA(LUN).EQ.0) GOTO 900
      IF(NTBB(LUN).EQ.0) GOTO 901
      IF(NTBD(LUN).EQ.0) GOTO 902

C  MAKE SURE EACH TABLE A ENTRY DEFINED AS A SEQUENCE
C  --------------------------------------------------

      DO I=1,NTBA(LUN)
      NEMO = TABA(I,LUN)(4:11)
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)
      IF(TAB.NE.'D') GOTO 903
      ENDDO

C  CHECK TABLE B CONTENTS
C  ----------------------

      DO ITAB=1,NTBB(LUN)
      CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)
      ENDDO

C  CHECK TABLE D CONTNETS
C  ----------------------

      DO ITAB=1,NTBD(LUN)
      CALL NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)
      ENDDO

      RETURN
900   CALL ABORT('CHEKSTAB - EMPTY TABLE A')
901   CALL ABORT('CHEKSTAB - EMPTY TABLE B')
902   CALL ABORT('CHEKSTAB - EMPTY TABLE D')
903   CALL ABORT('CHEKSTAB - NO SEQUENCE DEFINED FOR TABLE A: '//NEMO)
      END
