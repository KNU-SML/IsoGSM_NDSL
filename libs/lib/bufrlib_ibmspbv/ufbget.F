      SUBROUTINE UFBGET(LUNIT,TAB,I1,IRET,STR)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
 
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /USRBIT/ NBIT(INV_),MBIT(INV_)
      COMMON /SCALEU/ SCLU(-1000:1000)                                        
 
      CHARACTER*(*) STR
      CHARACTER*10  TAG,TGS(100)
      CHARACTER*8   SUBSET,BVAL
      CHARACTER*3   TYP
      DIMENSION     TAB(I1)
      EQUIVALENCE   (IVAL,RVAL),(JVAL,BVAL)
      LOGICAL       OPENIT
      REAL*8        VAL
      REAL*4 VALS
 
      DATA MAXTG /100/
 
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB,USRTPL,INVWIN)
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1
      UPS(NODE) = (IVAL+IRF(NODE))*SCLU(-ISC(NODE))
C-----------------------------------------------------------------------
 
      IRET = 0
      BVAL = ' '
 
C  MAKE SURE A FILE/MESSAGE IS OPEN FOR INPUT
C  ------------------------------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.GE.0) GOTO 900
      IF(IM.EQ.0) GOTO 901
 
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------
 
      IF(NSUB(LUN).EQ.MSUB(LUN)) THEN
         IRET = -1
         RETURN
      ENDIF
 
C  PARSE THE STRING
C  ----------------
 
      CALL STRING(STR,LUN,I1,0)
 
C  CREATE THE TEMPLATE FOR THE ZERO EXPANSION OF THIS SUBSET
C  ---------------------------------------------------------
 
      N = 1
      CALL USRTPL(LUN,N,N)
 
      NBIT(1) = 0
      MBIT(1) = MBYT(LUN)*8 + 16
 
      DO N=2,NVAL(LUN)
      NODE = INV(N,LUN)
      NBIT(N) = IBT(NODE)
      MBIT(N) = MBIT(N-1)+NBIT(N-1)
      IF(ITP(NODE).EQ.1 .OR. NODE.EQ.NODS(NNOD)) THEN
         NVAL(LUN) = N
         GOTO 20
      ENDIF
      ENDDO
20    CONTINUE
 
C  UNPACK ONLY THE NODES FOUND IN THE STRING
C  -----------------------------------------
 
      DO I=1,NNOD
      NODE = NODS(I)
      INVN = INVWIN(NODE,LUN,1,NVAL(LUN))
      IF(INVN.GT.0) THEN
         CALL UPBB(IVAL,NBIT(INVN),MBIT(INVN),LUN)
         IF(ITP(NODE).EQ.1) THEN
            TAB(I) = IVAL
         ELSEIF(ITP(NODE).EQ.2) THEN
            IF(IVAL.LT.MPS(NODE)) TAB(I) = UPS(NODE)
         ELSEIF(ITP(NODE).EQ.3) THEN
            IVAL = ISHFT(IVAL,64-NBIT(INVN))
            IVAL = IOR(IVAL,ISHFT(JVAL,-NBIT(INVN)))
            TAB(I) = RVAL
         ENDIF
      ELSE
         TAB(I) = 10E10
      ENDIF
      ENDDO
 
      RETURN
900   CALL ABT('FILE NOT OPEN FOR INPUT')
901   CALL ABT('NO MESSAGE OPEN        ')
      END
