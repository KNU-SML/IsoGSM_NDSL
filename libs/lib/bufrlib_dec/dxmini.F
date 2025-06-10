      SUBROUTINE DXMINI(LUN,MBAY,MBYT,MB4,MBA,MBB,MBD)
#include "ufbprm.h"
      integer TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),
     .                LD30(10),DXSTR(10)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      DIMENSION     MBAY(5000)

      DATA NMCODE/99/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  INITIALIZE THE MESSAGE
C  ----------------------

      MBIT = 0
      DO I=1,5000
      MBAY(I) = 0
      ENDDO

      IH   = 0
      ID   = 0
      IM   = 0
      IY   = 0
      MTYP = 11
      MSBT = IDXV
      NSUB = 1
      IDXS = MSBT+1
      LDXS = NXSTR(IDXS)

      NBY0 = 8
      NBY1 = 18
      NBY2 = 0
      NBY3 = 7 + NXSTR(IDXS) + 1
      NBY4 = 7
      NBY5 = 4
      MBYT = NBY0+NBY1+NBY2+NBY3+NBY4+NBY5

      IF(MOD(NBY3,2).NE.0) GOTO 900

C  SECTION 0
C  ---------

      CALL PKC('BUFR' ,  4 , MBAY,MBIT)
      CALL PKB(  MBYT , 24 , MBAY,MBIT)
      CALL PKB(     2 ,  8 , MBAY,MBIT)

C  SECTION 1
C  ---------

      CALL PKB(  NBY1 , 24 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)
      CALL PKB(NMCODE , 16 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)
      CALL PKB(  MTYP ,  8 , MBAY,MBIT)
      CALL PKB(  MSBT ,  8 , MBAY,MBIT)
      CALL PKB(     2 ,  8 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)
      CALL PKB(    IY ,  8 , MBAY,MBIT)
      CALL PKB(    IM ,  8 , MBAY,MBIT)
      CALL PKB(    ID ,  8 , MBAY,MBIT)
      CALL PKB(    IH ,  8 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)
      CALL PKB(     0 ,  8 , MBAY,MBIT)

C  SECTION 3
C  ---------

      CALL PKB(       NBY3 ,   24 , MBAY,MBIT)
      CALL PKB(          0 ,    8 , MBAY,MBIT)
      CALL PKB(          1 ,   16 , MBAY,MBIT)
      CALL PKB(       2**7 ,    8 , MBAY,MBIT)
C     CALL PKC(DXSTR(IDXS) , LDXS , MBAY,MBIT)
      do i=1,ldxs
      call pkb(iupm(dxstr(idxs)(i:i),8),8,mbay,mbit)
      enddo
      CALL PKB(          0 ,    8 , MBAY,MBIT)

C  SECTION 4
C  ---------

      MB4 = MBIT/8+1
      CALL PKB(NBY4 , 24 , MBAY,MBIT)
      CALL PKB(   0 ,  8 , MBAY,MBIT)
      MBA = MBIT/8+1
      CALL PKB(   0 ,  8 , MBAY,MBIT)
      MBB = MBIT/8+1
      CALL PKB(   0 ,  8 , MBAY,MBIT)
      MBD = MBIT/8+1
      CALL PKB(   0 ,  8 , MBAY,MBIT)

      IF(MBIT/8+NBY5.NE.MBYT) GOTO 901

      RETURN
900   CALL ABORT('DXMINI - UNEVEN SECTION 3')
901   CALL ABORT('DXMINI - BYTCNT IS OFF   ')
      END
