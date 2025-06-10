      SUBROUTINE UFBDMP(LUNIT,LUPRT)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
 
      CHARACTER*10 TAG,TG
      CHARACTER*8  VC
      CHARACTER*3  TYP,TP
      CHARACTER*1  YOU
      EQUIVALENCE  (VL,VC)
      REAL*8       VAL,VL
 
      DATA BMISS /10E10/
      DATA XMSG  /'MISSING'/
 
C----------------------------------------------------------------------
CFFP$ EXPAND (STATUS)
C----------------------------------------------------------------------
 
      if(luprt.eq.0) luout = 6
      if(luprt.ne.0) luout = luprt
 
C  CHECK THE FILE STATUS AND I-NODE
C  --------------------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IM.EQ.0) GOTO 901
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902
 
C  DUMP THE CONTENTS OF COMMON /USRINT/ FOR UNIT LUNIT
C  ---------------------------------------------------
 
      DO NV=1,NVAL(LUN)
      if(luprt.eq.0 .and. mod(nv,20).eq.0) then
         print*,'(MORE)'
         read(5,'(a1)') you
         if(you.eq.'q') return
      endif
      ND = INV (NV,LUN)
      VL = VAL (NV,LUN)
      TG = TAG (ND)
      TP = TYP (ND)
      IT = ITP (ND)
      IB = IBT (ND)
      IS = ISC (ND)
      IR = IRF (ND)
      JP = JUMP(ND)
      LK = LINK(ND)
      JB = JMPB(ND)
      RJ = RJUST(TG)
      IF(TP.NE.'CHR' .AND. NINT(VL).LT.BMISS) THEN
         WRITE(LUOUT,1) NV,TP,IT,TG,VL,IB,IS,IR,ND,JP,LK,JB
      ELSE
         IF(TP.NE.'CHR') VL = XMSG
         IF(TP.EQ.'CHR') VC = VC(1:IB/8)
         RJ = RJUST(VC)
         WRITE(LUOUT,2) NV,TP,IT,TG,VC,IB,IS,IR,ND,JP,LK,JB
      ENDIF
      ENDDO
 
1     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,F10.1,7(1X,I5))
2     FORMAT(I5,1X,A3,'-',I1,1X,A10,1X,A10  ,7(1X,I5))
 
 
C  EXITS
C  -----
 
      RETURN
900   CALL ABT('UFBDMP - FILE IS CLOSED                     ')
901   CALL ABT('UFBDMP - NO MESSAGE OPEN                    ')
902   CALL ABT('UFBDMP - I-NODE MISMATCH                    ')
      END
