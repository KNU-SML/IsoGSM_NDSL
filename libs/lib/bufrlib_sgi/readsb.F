      SUBROUTINE READSB(LUNIT,IRET)

      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)

C-----------------------------------------------------------------------
CFPP$ EXPAND(STATUS)
C-----------------------------------------------------------------------

      IRET = 0

C  CHECK THE FILE STATUS
C  ---------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) THEN
         IRET = -1
         RETURN
      ENDIF

C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------

      IF(NSUB(LUN).EQ.MSUB(LUN)) THEN
         IRET = -1
         RETURN
      ENDIF

C  READ THE NEXT SUBSET AND RESET THE POINTERS
C  -------------------------------------------

      IBIT = MBYT(LUN)*8
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)
      CALL RDTREE(LUN)
      MBYT(LUN) = MBYT(LUN) + NBYT
      NSUB(LUN) = NSUB(LUN) + 1

      RETURN
900   CALL ABORT('READSB - FILE IS CLOSED           ')
901   CALL ABORT('READSB - FILE IS OPEN FOR OUTPUT  ')
902   CALL ABORT('READSB - NO MESSAGE OPEN          ')
      END
