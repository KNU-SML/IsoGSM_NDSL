      SUBROUTINE PKTDD(ID,LUN,IDN,IRET)
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

C  ZERO THE COUNTER IF IDN IS ZERO
C  -------------------------------

      IF(IDN.EQ.0) THEN
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,0)
         IRET = 0
         RETURN
      ENDIF

C  UPDATE THE STORED DESCRIPTOR COUNT FOR THIS TABLE D ENTRY
C  ---------------------------------------------------------

      ND = IUPM(TABD(ID,LUN)(LDD:LDD),8)

      IF(ND.LT.0 .OR. ND.EQ.TBD_) THEN
         IRET = -1
         RETURN
      ELSE
         ND = ND+1
         CALL IPKM(TABD(ID,LUN)(LDD:LDD),1,ND)
         IRET = ND
      ENDIF

C  PACK AND STORE THE DESCRIPTOR
C  -----------------------------

      IDM = LDD+1 + (ND-1)*2
      CALL IPKM(TABD(ID,LUN)(IDM:IDM),2,IDN)

      RETURN
      END
