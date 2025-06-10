      SUBROUTINE UPTDD(ID,LUN,IENT,IRET)
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

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      LDD = LDXD(IDXV+1)+1

C  CHECK IF IENT IS IN BOUNDS
C  --------------------------

      NDSC = IUPM(TABD(ID,LUN)(LDD:LDD),8)

      IF(IENT.EQ.0) THEN
         IRET = NDSC
         RETURN
      ELSEIF(IENT.LT.0 .OR. IENT.GT.NDSC) THEN
         CALL ABT('UPTDD - IENT OUT OF RANGE')
      ENDIF

C  RETURN THE DESCRIPTOR INDICATED BY IENT
C  ---------------------------------------

      IDSC = LDD+1 + (IENT-1)*2
      IRET = IUPM(TABD(ID,LUN)(IDSC:IDSC),16)

      RETURN
      END
