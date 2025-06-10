      SUBROUTINE RDTREE(LUN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /USRBIT/ NBIT(INV_),MBIT(INV_)
      COMMON /SCALEU/ SCLU(-1000:1000)                                        

      CHARACTER*10 TAG
      CHARACTER*8  CVAL
      CHARACTER*3  TYP
      DIMENSION    IVAL(INV_)
      EQUIVALENCE  (CVAL,RVAL)
      REAL*8       VAL,RVAL,ups,ten
      data         ten/10./

C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1
c     UPS(NODE) = (IVAL(N)+IRF(NODE))*SCLU(-ISC(NODE))
      UPS(NODE) = (IVAL(N)+IRF(NODE))*ten**(-ISC(NODE))
C-----------------------------------------------------------------------

C  CYCLE THROUGH A SUBSET SETTING UP THE USER ARRAY
C  ------------------------------------------------

      MBIT(1) = IBIT
      NBIT(1) = 0
      CALL RCSTPL(LUN)

C  UNPACK A SUBSET INTO THE USER ARRAY
C  -----------------------------------

      DO N=1,NVAL(LUN)
      CALL UPBB(IVAL(N),NBIT(N),MBIT(N),LUN)
      ENDDO

C  CONVERT THE UNPACKED INTEGERS TO THE PROPER TYPES
C  -------------------------------------------------

      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.1) THEN
         VAL(N,LUN) = IVAL(N)
      ELSEIF(ITP(NODE).EQ.2) THEN
         IF(IVAL(N).LT.MPS(NODE)) VAL(N,LUN) = UPS(NODE)
      ENDIF
      ENDDO

C  SPECIAL TREATMENT FOR CHARACTERS
C  --------------------------------

      DO N=1,NVAL(LUN)
      NODE = INV(N,LUN)
      IF(ITP(NODE).EQ.3) THEN
         CVAL = ' '
         CALL UPC(CVAL,NBIT(N)/8,MBAY(1,LUN),MBIT(N))
         VAL(N,LUN) = RVAL
      ENDIF
      ENDDO

      RETURN
      END
