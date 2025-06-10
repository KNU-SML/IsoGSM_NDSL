      SUBROUTINE NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)
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
      CHARACTER*8   NEMO,NEMS,NEMT,NEMF
      CHARACTER*1   TAB
      DIMENSION     NEMS(TBD_),IRPS(TBD_),KNTS(TBD_)
      LOGICAL       REP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IF(ITAB.LE.0 .OR. ITAB.GT.NTBD(LUN)) GOTO 900

      REP  = .FALSE.

C  CLEAR THE RETURN VALUES
C  -----------------------

      NSEQ = 0

      DO I=1,TBD_
      NEMS(I) = ' '
      IRPS(I) = 0
      KNTS(I) = 0
      ENDDO

C  PARSE THE TABLE D ENTRY
C  -----------------------

      NEMO = TABD(ITAB,LUN)(7:14)
      IDSC = IDND(ITAB,LUN)
      CALL UPTDD(ITAB,LUN,0,NDSC)

      IF(IDSC.LT.IFXY('300000')) GOTO 901
      IF(IDSC.GT.IFXY('363255')) GOTO 901
C     IF(NDSC.LE.0             ) GOTO 902

      DO J=1,NDSC
      IF(NSEQ+1.GT.TBD_) GOTO 903
      CALL UPTDD(ITAB,LUN,J,IDSC)
      CALL NUMTAB(LUN,IDSC,NEMT,TAB,IRET)
      IF(TAB.EQ.'R') THEN
         IF(REP) GOTO 904
         REP = .TRUE.
         IF(IRET.LT.0) THEN
            IRPS(NSEQ+1) = 1
            KNTS(NSEQ+1) = ABS(IRET)
         ELSEIF(IRET.GT.0) THEN
            IRPS(NSEQ+1) = IRET
         ENDIF
      ELSEIF(TAB.EQ.'F') THEN
         IF(.NOT.REP) GOTO 904
         IRPS(NSEQ+1) = IRET
         REP = .FALSE.
      ELSEIF(TAB.EQ.'D') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         NEMS(NSEQ) = NEMT
      ELSEIF(TAB.EQ.'B') THEN
         REP = .FALSE.
         NSEQ = NSEQ+1
         IF(NEMT(1:1).EQ.'.') THEN
            CALL UPTDD(ITAB,LUN,J+1,IDSC)
            CALL NUMTAB(LUN,IDSC,NEMF,TAB,IRET)
            CALL RSVFVM(NEMT,NEMF)
            IF(TAB.NE.'B') GOTO 906
         ENDIF
         NEMS(NSEQ) = NEMT
      ELSE
         GOTO 905
      ENDIF
      ENDDO

      RETURN
900   CALL ABT('NEMTBD - ITAB NOT IN TABLE D   '                )
901   CALL ABT('NEMTBD - BAD DESCRIPTOR VALUE: '          //NEMO)
902   CALL ABT('NEMTBD - ZERO LENGTH SEQUENCE: '          //NEMO)
903   CALL ABT('NEMTBD - TOO MANY DESCRIPTORS IN SEQ: '   //NEMO)
904   CALL ABT('NEMTBD - REPLICATOR OUT OF ORDER IN SEQ: '//NEMO)
905   CALL ABT('NEMTBD - BAD DESCRIPTOR IN SEQUENCE: '    //NEMO)
906   CALL ABT('NEMTBD - FOLLOWING VALUE NOT FROM TABLEB:'//NEMF)
      END
