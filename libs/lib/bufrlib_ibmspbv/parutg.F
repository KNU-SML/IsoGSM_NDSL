      SUBROUTINE PARUTG(LUN,IO,UTG,NOD,KON,VAL,*)
#include "ufbprm.h"
      integer TAB_,SEQ_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /UTGPRM/ PICKEY

      CHARACTER*20  UTG,ATAG
      CHARACTER*10  TAG
      CHARACTER*3   TYP,ATYP,BTYP
      CHARACTER*1   COND(5)
      DIMENSION     BTYP(8),IOK(8)
      LOGICAL       PICKEY

      DATA NCHK   / 8/
      DATA BTYP   /'SUB','SEQ','REP','RPC','RPS','DRB','DRP','DRS'/
      DATA IOK    /  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,   0 ,   0 /
      DATA LTG    /20/

C----------------------------------------------------------------------
      PICKEY = .FALSE.
      COND(1) = '='
      COND(2) = '!'
      COND(3) = '<'
      COND(4) = '>'
      COND(5) = '^'
      NCOND   = 5
C----------------------------------------------------------------------

      ATAG  = ' '
      ATYP  = ' '
      KON   = 0
      NOD   = 0
      VAL   = 0

C  PARSE THE TAG
C  -------------

      DO I=1,LTG
      IF(UTG(I:I).EQ.' ') GOTO 1
      DO J=1,NCOND
      IF(UTG(I:I).EQ.COND(J)) THEN
         KON = J
         ICV = I+1
         GOTO 1
      ENDIF
      ENDDO
      ATAG(I:I) = UTG(I:I)
      ENDDO

C  FIND THE TAG IN THE SUBSET TABLE
C  --------------------------------

1     INOD = INODE(LUN)
      DO NOD=INOD,ISC(INOD)
      IF(ATAG.EQ.TAG(NOD)) GOTO 2
      ENDDO

      IF(KON.EQ.0 .AND. (IO.EQ.0.OR.ATAG.EQ.'NUL'.OR..NOT.PICKEY)) THEN
C     IF(KON.EQ.0) THEN
         NOD = 0
         RETURN
      ELSE
         PRINT*,'ABT TRYING TO WRITE A NON-EXISTANT MNEMONIC:'//ATAG
         RETURN 1
      ENDIF

C  CHECK FOR A VALID NODE TYP
C  --------------------------

2     IF(KON.EQ.5) THEN
         IF(TYP(NOD-1).NE.'DRP' .AND. TYP(NOD-1).NE.'DRS') GOTO 901
      ELSE
         ATYP = TYP(NOD)
         DO I=1,NCHK
         IF(ATYP.EQ.BTYP(I) .AND. IO.NE.IOK(I)) GOTO 902
         ENDDO
      ENDIF

C  IF A COND NODE GET THE COND VALUE
C  ---------------------------------

      IF(KON.NE.0) THEN
         CALL STRNUM(UTG(ICV:LTG),NUM)
         IF(NUM.LT.0) GOTO 903
         VAL = NUM
      ENDIF

      RETURN
900   CALL ABT('PARUTG - NO VALID TAG FOUND IN              :'//UTG)
901   CALL ABT('PARUTG - BUMP NODE MUST BE TYPE RPC(DRP)    :'//UTG)
902   CALL ABT('PARUTG - ILLEGAL NODE TYPE:'//ATYP//       ':'//UTG)
903   CALL ABT('PARUTG - BAD OR MISSING COND VALUE IN       :'//UTG)
      END
