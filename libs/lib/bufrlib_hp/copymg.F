      SUBROUTINE COPYMG(LUNIN,LUNOT)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*10 TAG
      CHARACTER*8  SUBSET
      CHARACTER*3  TYP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE FILE STATUSES
C  -----------------------

      CALL STATUS(LUNIN,LIN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.LE.0) GOTO 902

      CALL STATUS(LUNOT,LOT,IL,IM)
      IF(IL.EQ.0) GOTO 903
      IF(IL.LT.0) GOTO 904
      IF(IM.NE.0) GOTO 905

C  MAKE SURE BOTH FILES HAVE THE SAME TABLES
C  -----------------------------------------

      SUBSET = TAG(INODE(LIN))
      CALL NEMTBA(LOT,SUBSET,MSGT,MSTB,INOD)
      IF(INODE(LIN).NE.INOD) GOTO 906

C  EVERYTHING OKAY, COPY A MESSAGE
C  -------------------------------

      MBYM = IUPB(MBAY(1,LIN),5,24)
      CALL MSGWRT(LUNOT,MBAY(1,LIN),MBYM)

C  SET THE MESSAGE CONTROL WORLDS FOR LUNOT
C  ----------------------------------------

      NMSG (LOT) = NMSG(LOT) + 1
      NSUB (LOT) = MSUB(LIN)
      IDATE(LOT) = IDATE(LIN)
      INODE(LOT) = INODE(LIN)

C  NORMAL EXIT
C  -----------

      RETURN

C  ERROR EXITS
C  -----------

900   CALL ABORT('COPYMG - INPUT FILE IS CLOSED                    ')
901   CALL ABORT('COPYMG - INPUT FILE IS OPEN FOR OUTPUT           ')
902   CALL ABORT('COPYMG - NO INPUT FILE MESSAGE OPEN              ')
903   CALL ABORT('COPYMG - OUTPUT FILE IS CLOSED                   ')
904   CALL ABORT('COPYMG - OUTPUT FILE IS OPEN FOR OUTPUT          ')
905   CALL ABORT('COPYMG - OUTPUT FILE MESSAGE OPEN                ')
906   CALL ABORT('COPYMG - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES')
      END
