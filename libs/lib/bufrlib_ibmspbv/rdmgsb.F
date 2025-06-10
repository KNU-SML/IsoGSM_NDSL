      SUBROUTINE RDMGSB(LUNIT,IMSG,ISUB)

      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)

      CHARACTER*8 SUBSET

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  OPEN THE FILE AND SKIP TO MESSAGE IMSG
C  ----------------------------------------

      CALL OPENBF(LUNIT,'IN',LUNIT)
      CALL STATUS(LUNIT,LUN,IL,IM)

      DO I=1,IMSG-1
      READ(LUNIT,ERR=900,END=901)
      ENDDO

      CALL READMG(LUNIT,SUBSET,JDATE,IRET)
      IF(IRET.NE.0) GOTO 901

C  POSITION AT SUBSET ISUB
C  -------------------------

      DO I=1,ISUB-1
      IF(NSUB(LUN).GT.MSUB(LUN)) GOTO 902
      IBIT = MBYT(LUN)*8
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
      MBYT(LUN) = MBYT(LUN) + NBYT
      NSUB(LUN) = NSUB(LUN) + 1
      ENDDO

      CALL READSB(LUNIT,IRET)
      IF(IRET.NE.0) GOTO 902

C  NORMAL EXIT
C  -----------

      RETURN
900   CALL ABT('RDMGSB - ERROR READING MESSAGES  ')
901   CALL ABT('RDMGSB - EOF - READING MESSAGES  ')
902   CALL ABT('RDMGSB - EOM - READING SUBSETS   ')
      END
