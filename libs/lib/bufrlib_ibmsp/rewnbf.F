      SUBROUTINE REWNBF(LUNIT,ISR)

      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /BUFRSR/ JUNN,JILL,JIMM,JBIT,JBYT,JMSG,JSUB,KSUB,JNOD,JDAT,
     .                JSR(10),JBAY(5000)

      CHARACTER*20 MSTR
      CHARACTER*10 CUNIT

C     DATA JSR /10*0/

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      WRITE(CUNIT,'(A8,I2)') 'ON UNIT ',LUNIT

C  TRY TO TRAP BAD CALL PROBLEMS
C  -----------------------------

      IF(ISR.EQ.0) THEN
         CALL STATUS(LUNIT,LUN,IL,IM)
         IF(JSR(LUN).NE.0)  GOTO 900
         IF(IL.EQ.0) GOTO 901
      ELSEIF(ISR.EQ.1) THEN
         LUN = JUNN
         IF(JSR(JUNN).NE.1)  GOTO 902
      ELSE
         GOTO 903
      ENDIF

C  STORE FILE PARAMETERS AND SET FOR READING
C  -----------------------------------------

      IF(ISR.EQ.0) THEN
         JUNN = LUN
         JILL = IL
         JIMM = IM
         JBIT = IBIT
         JBYT = MBYT (LUN)
         JMSG = NMSG (LUN)
         JSUB = NSUB (LUN)
         KSUB = MSUB (LUN)
         JNOD = INODE(LUN)
         JDAT = IDATE(LUN)
         DO I=1,JBYT
         JBAY(I) = MBAY(I,LUN)
         ENDDO
         CALL WTSTAT(LUNIT,LUN,-1,0)
      ENDIF

C  REWIND THE FILE AND POSITION AFTER THE DICTIONARY
C  -------------------------------------------------

      REWIND LUNIT
1     READ(LUNIT,ERR=904,END=2) MSTR
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1
2     BACKSPACE LUNIT

C  RESTORE FILE PARAMETERS AND POSITION IT TO WHERE IT WAS SAVED
C  -------------------------------------------------------------

      IF(ISR.EQ.1) THEN
         LUN        = JUNN
         IL         = JILL
         IM         = JIMM
         IBIT       = JBIT
         MBYT (LUN) = JBYT
         NMSG (LUN) = JMSG
         NSUB (LUN) = JSUB
         MSUB (LUN) = KSUB
         INODE(LUN) = JNOD
         IDATE(LUN) = JDAT
         DO I=1,JBYT
         MBAY(I,LUN) = JBAY(I)
         ENDDO
         DO IMSG=1,JMSG
         READ(LUNIT,ERR=904,END=905)
         ENDDO
         CALL WTSTAT(LUNIT,LUN,IL,IM)
      ENDIF

C  NORMAL EXIT
C  -----------

      JSR(LUN) = MOD(JSR(LUN)+1,2)

      RETURN

C  ERROR ENCOUNTERED
C  -----------------

900   CALL ABT('REWNBF - ATTEMPT TO SAVE BEFORE RESTORE '//CUNIT)
901   CALL ABT('REWNBF - FILE IS CLOSED  '               //CUNIT)
902   CALL ABT('REWNBF - ATTEMPT TO RESTORE BEFORE SAVE '//CUNIT)
903   CALL ABT('REWNBF - SAVE/RESTORE FLAG NE 0 OR 1 '   //CUNIT)
904   CALL ABT('REWNBF - FILE READ ERROR '               //CUNIT)
905   CALL ABT('REWNBF - EOF TRYING TO REPOSTION FILE '  //CUNIT)
      END
