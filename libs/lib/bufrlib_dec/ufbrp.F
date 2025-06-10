      SUBROUTINE UFBRP(LUN,USR,I1,I2,IO,IRET)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
 
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(TAB_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10)
 
      CHARACTER*10 TAG
      CHARACTER*3  TYP
      REAL*8       USR(I1,I2),VAL
      REAL*4 VALS
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IRET = 0
      INS1 = 0
      INS2 = 0
 
C  FRAME A SECTION OF THE BUFFER - RETURN WHEN NO FRAME
C  ----------------------------------------------------
 
1     IF(INS1+1.GT.NVAL(LUN)) RETURN
      INS1 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))
      IF(INS1.EQ.0) RETURN
 
      INS2 = INVtag(NODS(1),LUN,INS1+1,NVAL(LUN))
      IF(INS2.EQ.0) INS2 = NVAL(LUN)
      IRET = IRET+1
 
C  READ USER VALUES
C  ----------------
 
      IF(IO.EQ.0 .AND. IRET.LE.I2) THEN
         DO I=1,NNOD
         USR(I,IRET) = 10E10
         IF(NODS(I).GT.0) THEN
            INVN = INVtag(NODS(I),LUN,INS1,INS2)
            IF(INVN.GT.0) USR(I,IRET) = VAL(INVN,LUN)
         ENDIF
         ENDDO
      ENDIF
 
C  WRITE USER VALUES
C  -----------------
 
      IF(IO.EQ.1 .AND. IRET.LE.I2) THEN
         DO I=1,NNOD
         IF(NODS(I).GT.0) THEN
            INVN = INVtag(NODS(I),LUN,INS1,INS2)
            IF(INVN.GT.0) VAL(INVN,LUN) = USR(I,IRET)
         ENDIF
         ENDDO
      ENDIF
 
C  GO FOR NEXT FRAME
C  -----------------
 
      GOTO 1
 
      END
