      SUBROUTINE WRITDX(LUNIT,LUN,LUNDX)
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
      CHARACTER*6   ADN30
      CHARACTER*1   MOCT(24000)
      DIMENSION     MBAY(5000)
      EQUIVALENCE   (MOCT(1),MBAY(1))

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK UNITS, GET A DX TABLE, AND START A DX MESSAGE
C  ---------------------------------------------------

      IF(LUNIT.EQ.LUNDX) GOTO 900
      CALL READDX(LUNIT,LUN,LUNDX)
      CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)

      LDA = LDXA(IDXV+1)
      LDB = LDXB(IDXV+1)
      LDD = LDXD(IDXV+1)
      L30 = LD30(IDXV+1)

C  COPY TABLE A CONTENTS TO A BUFR DX MESSAGE
C  ------------------------------------------

      DO I=1,NTBA(LUN)
      IF(MBYT+LDA.GT.MAXDX) THEN
         CALL MSGWRT(LUNIT,MBAY,MBYT)
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)
      ENDIF
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LDA)
      CALL IPKM(MOCT(MBYA),1,IUPM(MOCT(MBYA), 8)+  1)
      MBIT = 8*(MBYB-1)
      CALL PKC(TABA(I,LUN),LDA,MBAY,MBIT)
      CALL PKB(          0,  8,MBAY,MBIT)
      CALL PKB(          0,  8,MBAY,MBIT)
      MBYT = MBYT+LDA
      MBYB = MBYB+LDA
      MBYD = MBYD+LDA
      ENDDO

C  COPY TABLE B CONTENTS TO A BUFR DX MESSAGE
C  ------------------------------------------

      DO I=1,NTBB(LUN)
      IF(MBYT+LDB.GT.MAXDX) THEN
         CALL MSGWRT(LUNIT,MBAY,MBYT)
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)
      ENDIF
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LDB)
      CALL IPKM(MOCT(MBYB),1,IUPM(MOCT(MBYB), 8)+  1)
      MBIT = 8*(MBYD-1)
      CALL PKC(TABB(I,LUN),LDB,MBAY,MBIT)
      CALL PKB(          0,  8,MBAY,MBIT)
      MBYT = MBYT+LDB
      MBYD = MBYD+LDB
      ENDDO

C  COPY TABLE D CONTENTS TO A BUFR DX MESSAGE
C  ------------------------------------------

      DO I=1,NTBD(LUN)
      NSEQ = IUPM(TABD(I,LUN)(LDD+1:LDD+1),8)
      LEND = LDD+1 + L30*NSEQ
      IF(MBYT+LEND.GT.MAXDX) THEN
         CALL MSGWRT(LUNIT,MBAY,MBYT)
         CALL DXMINI(LUN,MBAY,MBYT,MBY4,MBYA,MBYB,MBYD)
      ENDIF
      CALL IPKM(MOCT(MBY4),3,IUPM(MOCT(MBY4),24)+LEND)
      CALL IPKM(MOCT(MBYD),1,IUPM(MOCT(MBYD), 8)+   1)
      MBIT = 8*(MBYT-4)
      CALL PKC(TABD(I,LUN),LDD,MBAY,MBIT)
      CALL PKB(       NSEQ,  8,MBAY,MBIT)
         DO J=1,NSEQ
         JJ  = LDD+2 + (J-1)*2
         IDN = IUPM(TABD(I,LUN)(JJ:JJ),16)
         CALL PKC(ADN30(IDN,L30),L30,MBAY,MBIT)
         ENDDO
      MBYT = MBYT+LEND
      ENDDO

C  WRITE THE UNWRITTEN MESSAGE
C  ---------------------------

      CALL MSGWRT(LUNIT,MBAY,MBYT)

      RETURN
900   CALL ABORT('WRITDX - INPUR AND OUTPUT UNIT MUST NOT BE THE SAME')
      END
