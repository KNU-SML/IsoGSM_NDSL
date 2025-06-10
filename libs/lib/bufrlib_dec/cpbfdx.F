      SUBROUTINE CPBFDX(LUD,LUN)
#include "ufbprm.h"
      integer TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  INITIALIZE THE DX-TABLE PARTITION
C  ---------------------------------

      CALL DXINIT(LUN,0)

C  COPY ONE TABLE PARTITION TO ANOTHER
C  -----------------------------------

      INODE(LUN) = INODE(LUD)

      NTBA(LUN) = NTBA(LUD)
      NTBB(LUN) = NTBB(LUD)
      NTBD(LUN) = NTBD(LUD)

      DO I=1,NTBA(LUD)
      IDNA(I,LUN,1) = IDNA(I,LUD,1)
      IDNA(I,LUN,2) = IDNA(I,LUD,2)
      TABA(I,LUN) = TABA(I,LUD)
      MTAB(I,LUN) = MTAB(I,LUD)
      ENDDO

      DO I=1,NTBB(LUD)
      IDNB(I,LUN) = IDNB(I,LUD)
      TABB(I,LUN) = TABB(I,LUD)
      ENDDO

      DO I=1,NTBD(LUD)
      IDND(I,LUN) = IDND(I,LUD)
      TABD(I,LUN) = TABD(I,LUD)
      ENDDO

      RETURN
      END
