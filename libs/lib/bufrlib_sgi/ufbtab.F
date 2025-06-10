      SUBROUTINE UFBTAB(LUNIT,TAB,I1,I2,IRET,STR)
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
      COMMON /SCALEU/ SCLU(-1000:1000)                                        

      CHARACTER*(*) STR
      CHARACTER*10  TAG,TGS(100)
      CHARACTER*8   SUBSET,CVAL
      CHARACTER*3   TYP
      DIMENSION     TAB(I1,I2),JBAY(5000)
      EQUIVALENCE   (CVAL,RVAL)
      LOGICAL       OPENIT
      REAL*8        VAL,TAB
      REAL*4 VALS

      DATA MAXTG /100/

C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB,USRTPL)
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1
      UPS(NODE) = (IVAL+IRF(NODE))*SCLU(-ISC(NODE))
C-----------------------------------------------------------------------

      IRET = 0
      IREC = 0
      ISUB = 0

      DO J=1,I2
      DO I=1,I1
      TAB(I,J) = 10E10
      ENDDO
      ENDDO

C  SEE IF WE NEED TO OPEN A FILE
C  -----------------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      OPENIT = IL.EQ.0

      IF(OPENIT) THEN
         CALL OPENBF(LUNIT,'IN',LUNIT)
      ELSE
         CALL REWNBF(LUNIT,0)
      ENDIF

C  CHECK FOR SPECIAL TAGS IN STRING
C  --------------------------------

      CALL PARSEQ(STR,TGS,MAXTG,NTG)
      DO I=1,NTG
      IF(TGS(I).EQ.'IREC') IREC = I
      IF(TGS(I).EQ.'ISUB') ISUB = I
      ENDDO

C  READ A MESSAGE AND PARSE A STRING
C  ---------------------------------

10    CALL READMG(LUNIT,SUBSET,IDATE,MRET)
      IF(MRET.NE.0) GOTO 25
      CALL STRING(STR,LUN,I1,0)
      IF(IREC.GT.0) NODS(IREC) = 0
      IF(ISUB.GT.0) NODS(ISUB) = 0

C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------

15    IF(NSUB(LUN).EQ.MSUB(LUN)) GOTO 10
      IF(IRET+1.GT.I2) CALL ABORT('UFBTAB - TAB TOO SMALL')
      IRET = IRET+1

      DO I=1,NNOD
      NODS(I) = ABS(NODS(I))
      ENDDO

C  PARSE THE STRING NODES FROM A SUBSET
C  ------------------------------------

      MBIT = MBYT(LUN)*8 + 16
      NBIT = 0
      N = 1
      CALL USRTPL(LUN,N,N)
20    IF(N+1.LE.NVAL(LUN)) THEN
         N = N+1
         NODE = INV(N,LUN)
         MBIT = MBIT+NBIT
         NBIT = IBT(NODE)
         IF(ITP(NODE).EQ.1) THEN
            CALL UPBB(IVAL,NBIT,MBIT,LUN)
            CALL USRTPL(LUN,N,IVAL)
         ENDIF
         DO I=1,NNOD
         IF(NODS(I).EQ.NODE) THEN
            CALL UPBB(IVAL,NBIT,MBIT,LUN)
            IF(ITP(NODE).EQ.1) THEN
               TAB(I,IRET) = IVAL
            ELSEIF(ITP(NODE).EQ.2) THEN
               IF(IVAL.LT.MPS(NODE)) TAB(I,IRET) = UPS(NODE)
            ELSEIF(ITP(NODE).EQ.3) THEN
               CVAL = ' '
               CALL UPC(CVAL,NBIT/8,MBAY(1,LUN),MBIT)
               TAB(I,IRET) = RVAL
            ENDIF
            NODS(I) = -NODS(I)
            GOTO 20
         ENDIF
         ENDDO
         DO I=1,NNOD
         IF(NODS(I).GT.0) GOTO 20
         ENDDO
      ENDIF

C  UPDATE THE SUBSET POINTERS BEFORE NEXT READ
C  -------------------------------------------

      IBIT = MBYT(LUN)*8
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
      MBYT(LUN) = MBYT(LUN) + NBYT
      NSUB(LUN) = NSUB(LUN) + 1
      IF(IREC.GT.0) TAB(IREC,IRET) = NMSG(LUN)
      IF(ISUB.GT.0) TAB(ISUB,IRET) = NSUB(LUN)
      GOTO 15

C  LEAVE THE FILE AS IT WAS BEFORE
C  -------------------------------

25    IF(OPENIT) THEN
         CALL CLOSBF(LUNIT)
      ELSE
         CALL REWNBF(LUNIT,1)
      ENDIF

      RETURN
      END
