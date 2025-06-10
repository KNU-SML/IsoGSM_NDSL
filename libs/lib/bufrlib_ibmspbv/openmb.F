      SUBROUTINE OPENMB(LUNIT,SUBSET,JDATE)

      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)

      CHARACTER*(*) SUBSET
      LOGICAL       OPEN

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE FILE STATUS
C  ---------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.LT.0) GOTO 901

C  GET SOME SUBSET PARTICULARS
C  ---------------------------

      CALL NEMTBA(LUN,SUBSET,MTYP,MSTB,INOD)
      OPEN = IM.EQ.0 .OR. INOD.NE.INODE(LUN) .OR. JDATE.NE.IDATE(LUN)

C  MAYBE OPEN A NEW OR DIFFERENT TYPE OF MESSAGE
C  ---------------------------------------------

      IF(OPEN) THEN
         CALL CLOSMG(LUNIT)
         CALL WTSTAT(LUNIT,LUN,IL, 1)
         INODE(LUN) = INOD
         IDATE(LUN) = JDATE
         CALL MSGINI(LUN)
         CALL USRTPL(LUN,1,1)
         CALL STRCLN
      ENDIF

      RETURN
900   CALL ABT('OPENMB - FILE IS CLOSED            ')
901   CALL ABT('OPENMB - FILE IS OPEN FOR INPUT    ')
      END
