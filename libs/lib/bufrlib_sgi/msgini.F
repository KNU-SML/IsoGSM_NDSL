      SUBROUTINE MSGINI(LUN)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)

      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4
      COMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*10 TAG
      CHARACTER*8  SUBTAG
      CHARACTER*4  BUFR,SEVN
      CHARACTER*3  TYP
      CHARACTER*1  TAB

      DATA BUFR/'BUFR'/
      DATA SEVN/'7777'/
      DATA NMCD/   99 /

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  GET THE MESSAGE TAG AND TYPE, AND BREAK UP THE DATE
C  ---------------------------------------------------

      SUBTAG = TAG(INODE(LUN))
      CALL NEMTBA(LUN,SUBTAG,MTYP,MSBT,INOD)
      IF(INODE(LUN).NE.INOD) GOTO 900
      CALL NEMTAB(LUN,SUBTAG,ISUB,TAB,IRET)
      IF(IRET.EQ.0) GOTO 901

      MEAR = MOD(IDATE(LUN)/1000000,100)
      MMON = MOD(IDATE(LUN)/10000  ,100)
      MDAY = MOD(IDATE(LUN)/100    ,100)
      MOUR = MOD(IDATE(LUN)/1      ,100)

C  INITIALIZE THE MESSAGE
C  ----------------------

      MBIT = 0
      NBY0 = 8
      NBY1 = 18
      NBY2 = 0
      NBY3 = 20
      NBY4 = 4
      NBY5 = 4
      NBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5

C  SECTION 0
C  ---------

      CALL PKC(BUFR ,  4 , MBAY(1,LUN),MBIT)
      CALL PKB(NBYT , 24 , MBAY(1,LUN),MBIT)
      CALL PKB(   2 ,  8 , MBAY(1,LUN),MBIT)

C  SECTION 1
C  ---------

      CALL PKB(NBY1 , 24 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(NMCD , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MTYP ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MSBT ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   2 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MEAR ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MMON ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MDAY ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(MOUR ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)

C  SECTION 3
C  ---------

      CALL PKB(NBY3 , 24 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(2**7 ,  8 , MBAY(1,LUN),MBIT)
      CALL PKB(IBCT , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(ISUB , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(IPD1 , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(IPD2 , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(IPD3 , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(IPD4 , 16 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)

C  SECTION 4
C  ---------

      CALL PKB(NBY4 , 24 , MBAY(1,LUN),MBIT)
      CALL PKB(   0 ,  8 , MBAY(1,LUN),MBIT)

C  SECTION 5
C  ---------

      CALL PKC(SEVN ,  4 , MBAY(1,LUN),MBIT)

C  DOUBLE CHECK INITIAL MESSAGE LENGTH
C  -----------------------------------

      IF(MOD(MBIT,8).NE.0) GOTO 902
      IF(MBIT/8.NE.NBYT  ) GOTO 903

      NMSG(LUN) = NMSG(LUN)+1
      NSUB(LUN) = 0
      MBYT(LUN) = NBYT

      RETURN
900   CALL ABORT('MSGINI - INODE MISMATCH                      ')
901   CALL ABORT('MSGINI - CANT FIND SUBSET IDN: '//SUBTAG      )
902   CALL ABORT('MSGINI - INITIAL MESSAGE OUT OF BYTE BOUNDARY')
903   CALL ABORT('MSGINI - INITIAL MESSAGE COUNT FAILS CHECK   ')
      END
