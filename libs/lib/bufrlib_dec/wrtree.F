      SUBROUTINE WRTREE(LUN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /SCALEU/ SCLU(-1000:1000)                                        

      CHARACTER*10 TAG
      CHARACTER*8  CVAL
      CHARACTER*3  TYP
      DIMENSION    IVAL(INV_)
      EQUIVALENCE  (CVAL,RVAL)
      REAL*8       VAL,RVAL,pks,ten
      data ten /10./

C-----------------------------------------------------------------------
c     PKS(NODE) = VAL(N,LUN)*SCLU(ISC(NODE))-IRF(NODE)
      PKS(NODE) = VAL(N,LUN)*ten**(ISC(NODE))-IRF(NODE)
C-----------------------------------------------------------------------

C  CONVERT USER NUMBERS INTO SCALED INTEGERS
C  -----------------------------------------

      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.1) THEN
         IVAL(N) = VAL(N,LUN)
      ELSEIF(ITP(NODE).EQ.2) THEN
         IVAL(N) = NINT(PKS(NODE))
      ENDIF
      ENDDO

C  PACK THE USER ARRAY INTO THE SUBSET BUFFER
C  ------------------------------------------

      IBIT = 16

      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).LT.3) THEN
         CALL PKB(IVAL(N),IBT(NODE),IBAY,IBIT)
      ELSE
         RVAL = VAL(N,LUN)
         CALL PKC(CVAL,IBT(NODE)/8,IBAY,IBIT)
      ENDIF
      ENDDO

      RETURN
      END
