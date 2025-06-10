      SUBROUTINE UFBRW(LUN,USR,I1,I2,IO,IRET)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10)

      CHARACTER*10 TAG
      CHARACTER*3  TYP
      REAL*8       USR(I1,I2),VAL
      REAL*4 VALS

C----------------------------------------------------------------------
CFPP$ EXPAND (CONWIN,DRSTPL,GETWIN,INVWIN,LSTRPS,NEWIN,NXTWIN)
C----------------------------------------------------------------------

      IRET = 0

C  LOOP OVER COND WINDOWS
C  ----------------------

      INC1 = 1
      INC2 = 1

1     CALL CONWIN(LUN,INC1,INC2,I2)
      IF(NNOD.EQ.0) THEN
         IRET = I2
         RETURN
      ELSEIF(INC1.EQ.0) THEN
         RETURN
      ELSE
         DO I=1,NNOD
         IF(NODS(I).GT.0) THEN
            INS2 = INC1
            CALL GETWIN(NODS(I),LUN,INS1,INS2)
            IF(INS1.EQ.0) RETURN
            GOTO 2
         ENDIF
         ENDDO
         INS1 = INC1
         INS2 = INC2
      ENDIF

C  LOOP OVER STORE NODES
C  ---------------------

2     IRET = IRET+1

C     print*,'ufbrw:',iret,':',ins1,':',ins2,':',inc1,':',inc2
C     print'(5a10)',(tag(inv(i,lun)),i=ins1,ins2)

C  WRITE USER VALUES
C  -----------------

      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN
         DO I=1,NNOD
         IF(NODS(I).GT.0) THEN
            IF(USR(I,IRET).NE.10E10) THEN
               INVN = INVWIN(NODS(I),LUN,INS1,INS2)
               IF(INVN.EQ.0) THEN
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)
                  if(invn.eq.0) then
                     iret = 0
                     return
                  endif
                  CALL NEWWIN(LUN,INC1,INC2)
                  VAL(INVN,LUN) = USR(I,IRET)
               ELSEIF(LSTRPS(NODS(I),LUN).EQ.0) THEN
                  VAL(INVN,LUN) = USR(I,IRET)
               ELSEIF(VAL(INVN,LUN).EQ.10E10) THEN
                  VAL(INVN,LUN) = USR(I,IRET)
               ELSE
                  CALL DRSTPL(NODS(I),LUN,INS1,INS2,INVN)
                  if(invn.eq.0) then
                     iret = 0
                     return
                  endif
                  CALL NEWWIN(LUN,INC1,INC2)
                  VAL(INVN,LUN) = USR(I,IRET)
               ENDIF
            ENDIF
         ENDIF
         ENDDO
      ENDIF

C  READ USER VALUES
C  ----------------

      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN
         DO I=1,NNOD
         USR(I,IRET) = 10E10
         IF(NODS(I).GT.0) THEN
            INVN = INVWIN(NODS(I),LUN,INS1,INS2)
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)
         ENDIF
         ENDDO
      ENDIF

C  DECIDE WHAT TO DO NEXT
C  ----------------------

      IF(IO.EQ.1.AND.IRET.EQ.I2) RETURN
      CALL NXTWIN(LUN,INS1,INS2)
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2
      IF(NCON.GT.0) GOTO 1

      RETURN
      END
