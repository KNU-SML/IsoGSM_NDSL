      SUBROUTINE MAKESTAB
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)

      COMMON /QUIET/  IPRT
      COMMON /STBFR/  IOLUN(10),IOMSG(10)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      CHARACTER*10  TAG
      CHARACTER*8   NEMO
      CHARACTER*3   TYP
      DIMENSION     LUS(10)
      LOGICAL       EXPAND,PRTTAB
      REAL*8        VAL

      DATA PRTTAB /.FALSE./

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      PRTTAB = IPRT.GE.2

C  RESET POINTER TABLE AND STRING CACHE
C  ------------------------------------

      NTAB = 0
      CALL STRCLN

C  FIGURE OUT WHICH UNITS SHARE TABLES
C  -----------------------------------

      DO LUN=1,10
      LUS(LUN) = 0
      IF(IOLUN(LUN).NE.0) THEN
         DO LUM=1,LUN-1
         IF(MTAB(1,LUN).EQ.MTAB(1,LUM)) LUS(LUN) = LUM
         ENDDO
      ENDIF
      ENDDO

C  INITIALIZE JUMP-LINK TABLES WITH SUBSETS/SEQUENCES/ELEMENTS
C  ----------------------------------------------------------

      DO LUN=1,10

      IF(IOLUN(LUN).NE.0) THEN

C  RESET ANY EXISTING INVENTORY POINTERS
C  -------------------------------------

         IF(IOMSG(LUN).NE.0) THEN
            IF(LUS(LUN).EQ.0) INC = (NTAB+1)-MTAB(1,LUN)
            IF(LUS(LUN).NE.0) INC = MTAB(1,LUS(LUN))-MTAB(1,LUN)
            DO N=1,NVAL(LUN)
            INV(N,LUN) = INV(N,LUN)+INC
            ENDDO
         ENDIF

C  CREATE NEW TABLE ENTRIES IF THIS UNIT DOESNT SHARE EXISTING ONES
C  -----------------------------------------------------------------

         IF(LUS(LUN).EQ.0) THEN
            CALL CHEKSTAB(LUN)
            DO ITBA=1,NTBA(LUN)
            INOD = NTAB+1
            NEMO = TABA(ITBA,LUN)(4:11)
            CALL TABSUB(LUN,NEMO)
            MTAB(ITBA,LUN) = INOD
            ISC(INOD)      = NTAB
C           DO N1=INOD,ISC(INOD)-1
C           DO N2=N1+1,ISC(INOD)
C           IF(TAG(N1).EQ.TAG(N2)) GOTO 900
C           ENDDO
C           ENDDO
            ENDDO
         ENDIF

      ENDIF
      ENDDO

C  STORE TYPES AND INITIAL VALUES AND COUNTS
C  -----------------------------------------

      DO NODE=1,NTAB
      IF(TYP(NODE).EQ.'SUB') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 1
         ITP (NODE) = 0
      ELSEIF(TYP(NODE).EQ.'SEQ') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 1
         ITP (NODE) = 0
      ELSEIF(TYP(NODE).EQ.'RPC') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 0
         ITP (NODE) = 0
      ELSEIF(TYP(NODE).EQ.'RPS') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 0
         ITP (NODE) = 0
      ELSEIF(TYP(NODE).EQ.'REP') THEN
         VALI(NODE) = 10E10
         KNTI(NODE) = IRF(NODE)
         ITP (NODE) = 0
      ELSEIF(TYP(NODE).EQ.'DRS') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 1
         ITP (NODE) = 1
      ELSEIF(TYP(NODE).EQ.'DRP') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 1
         ITP (NODE) = 1
      ELSEIF(TYP(NODE).EQ.'DRB') THEN
         VALI(NODE) = 0
         KNTI(NODE) = 0
         ITP (NODE) = 1
      ELSEIF(TYP(NODE).EQ.'NUM') THEN
         VALI(NODE) = 10E10
         KNTI(NODE) = 1
         ITP (NODE) = 2
      ELSEIF(TYP(NODE).EQ.'CHR') THEN
         VALI(NODE) = 10E10
         KNTI(NODE) = 1
         ITP (NODE) = 3
      ELSE
         GOTO 901
      ENDIF
      ENDDO

C  SET UP EXPANSION SEGMENTS FOR TYPE 'SUB', 'DRP', AND 'DRS' NODES
C  ----------------------------------------------------------------

      NEWN = 0

      DO N=1,NTAB
      ISEQ(N,1) = 0
      ISEQ(N,2) = 0
      EXPAND = TYP(N).EQ.'SUB' .OR. TYP(N).EQ.'DRP' .OR. TYP(N).EQ.'DRS'
     .                         .OR. TYP(N).EQ.'REP' .OR. TYP(N).EQ.'DRB'
      IF(EXPAND) THEN
         ISEQ(N,1) = NEWN+1
         NODA = N
         NODE = N+1
         DO K=1,TAB_
         KNT(K) = 0
         ENDDO
         IF(TYP(NODA).EQ.'REP') KNT(NODE) = KNTI(NODA)
         IF(TYP(NODA).NE.'REP') KNT(NODE) = 1

1        NEWN = NEWN+1
         IF(NEWN.GT.SEQ_) GOTO 902
         JSEQ(NEWN) = NODE
         KNT(NODE) = MAX(KNTI(NODE),KNT(NODE))
2        IF(JUMP(NODE)*KNT(NODE).GT.0) THEN
            NODE = JUMP(NODE)
            GOTO 1
         ELSE IF(LINK(NODE).GT.0) THEN
            NODE = LINK(NODE)
            GOTO 1
         ELSE
            NODE = JMPB(NODE)
            IF(NODE.EQ.NODA) GOTO 3
            IF(NODE.EQ.0   ) GOTO 903
            KNT(NODE) = MAX(KNT(NODE)-1,0)
            GOTO 2
         ENDIF
3        ISEQ(N,2) = NEWN
      ENDIF
      ENDDO

C  PRINT THE SEQUENCE TABLES
C  ------------------------

      IF(PRTTAB) THEN
         PRINT*
         DO I=1,NTAB
         PRINT99,I,
     .   TAG(I),TYP(I),JMPB(I),JUMP(I),LINK(I),IBT(I),IRF(I),ISC(I)
         ENDDO
         PRINT*
99       FORMAT(I5,2X,A10,A5,6I8)
      ENDIF

      RETURN
900   CALL ABORT('MAKESTAB - DUP IN SUBSET: '//TAG(N1)//':'//NEMO)
901   CALL ABORT('MAKESTAB - UNKNOWN TYPE : '         //TYP(NODE))
902   CALL ABORT('MAKESTAB - JSEQ OVERFLOW       : '  //TAG(N   ))
903   CALL ABORT('MAKESTAB - FAILED TO CIRCULATE : '  //TAG(N   ))
      END
