      SUBROUTINE NEMTAB(LUN,NEMO,IDN,TAB,IRET)
#include "ufbprm.h"
      integer TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
 
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)
 
      CHARACTER*(*) NEMO
      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      CHARACTER*8   NEMT
      CHARACTER*1   TAB
      LOGICAL       FOLVAL
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      FOLVAL = NEMO(1:1).EQ.'.'
      IRET = 0
      TAB = ' '
 
C  LOOK FOR NEMO IN TABLE B
C  ------------------------
 
      DO 1 I=1,NTBB(LUN)
      NEMT = TABB(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ELSEIF(FOLVAL.AND.NEMT(1:1).EQ.'.') THEN
         DO J=2,LEN(NEMT)
         IF(NEMT(J:J).NE.'.' .AND. NEMT(J:J).NE.NEMO(J:J)) GOTO 1
         ENDDO
         IDN  = IDNB(I,LUN)
         TAB  = 'B'
         IRET = I
         RETURN
      ENDIF
1     ENDDO

C  DONT LOOK IN TABLE D FOR FOLLOWING VALUE-MNEMONICS
C  ---------------------------------------------------

      IF(FOLVAL) RETURN
 
C  LOOK IN TABLE D IF WE GOT THIS FAR
C  ----------------------------------
 
      DO I=1,NTBD(LUN)
      NEMT = TABD(I,LUN)(7:14)
      IF(NEMT.EQ.NEMO) THEN
         IDN  = IDND(I,LUN)
         TAB  = 'D'
         IRET = I
         RETURN
      ENDIF
      ENDDO
 
      RETURN
      END
