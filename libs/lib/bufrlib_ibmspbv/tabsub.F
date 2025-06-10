      SUBROUTINE TABSUB(LUN,NEMO)
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)

      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*10 TAG
      CHARACTER*8  NEMO,NEMS,NEM
      CHARACTER*3  TYP
      CHARACTER*1  TAB
      DIMENSION    NEM(TBD_,10),IRP(TBD_,10),KRP(TBD_,10)
      DIMENSION    DROP(10),JMP0(10),NODL(10),NTAG(10,2)
      LOGICAL      DROP

      DATA MAXLIM /10/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MNEMONIC
C  ------------------

      CALL NEMTAB(LUN,NEMO,IDN,TAB,ITAB)
      IF(TAB.NE.'D') GOTO 900

C  STORE A SUBSET NODE AND JUMP/LINK THE TREE
C  ------------------------------------------

      CALL INCTAB(NEMO,'SUB',NODE)
      JUMP(NODE) = NODE+1
      JMPB(NODE) = 0
      LINK(NODE) = 0
      IBT (NODE) = 0
      IRF (NODE) = 0
      ISC (NODE) = 0

      CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,1),IRP(1,1),KRP(1,1))
      NTAG(1,1) = 1
      NTAG(1,2) = NSEQ
      JMP0(1)   = NODE
      LIMB      = 1

C  THIS LOOP RESOLVES ENTITIES IN A SUBSET BY EMULATING RECURSION
C  --------------------------------------------------------------

1     DO N=NTAG(LIMB,1),NTAG(LIMB,2)

      NTAG(LIMB,1) = N+1
      NODL(LIMB)   = NTAB+1
      DROP(LIMB)   = N.EQ.NTAG(LIMB,2)

      CALL NEMTAB(LUN,NEM(N,LIMB),IDN,TAB,ITAB)
      NEMS = NEM(N,LIMB)
      IREP = IRP(N,LIMB)
      IKNT = KRP(N,LIMB)
      JUM0 = JMP0(LIMB)
      CALL TABENT(LUN,NEMS,TAB,ITAB,IREP,IKNT,JUM0)

      IF(TAB.EQ.'D') THEN
         LIMB = LIMB+1
         IF(LIMB.GT.MAXLIM) GOTO 901
         CALL NEMTBD(LUN,ITAB,NSEQ,NEM(1,LIMB),IRP(1,LIMB),KRP(1,LIMB))
         NTAG(LIMB,1) = 1
         NTAG(LIMB,2) = NSEQ
         JMP0(LIMB)   = NTAB
         GOTO 1
      ELSE IF(DROP(LIMB)) THEN
2        LINK(NODL(LIMB)) = 0
         LIMB = LIMB-1
         IF(LIMB.EQ.0 ) RETURN
         IF(DROP(LIMB)) GOTO 2
         LINK(NODL(LIMB)) = NTAB+1
         GOTO 1
      ELSE
         LINK(NODL(LIMB)) = NTAB+1
      ENDIF

      ENDDO

      CALL ABT('TABSUB - SHOULD NOT GET HERE               ')
900   CALL ABT('TABSUB - SUBSET NODE NOT IN TABLE D: '//NEMO)
901   CALL ABT('TABSUB - TOO MANY LIMBS                    ')
      END
