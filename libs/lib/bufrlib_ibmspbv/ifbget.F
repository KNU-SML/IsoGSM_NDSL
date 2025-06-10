      FUNCTION IFBGET(LUNIT)
 
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  MAKE SURE A FILE/MESSAGE IS OPEN FOR INPUT
C  ------------------------------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.GE.0) GOTO 900
      IF(IM.EQ.0) GOTO 901
 
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE
C  ---------------------------------------------
 
      IF(NSUB(LUN).LT.MSUB(LUN)) THEN
         IFBGET = 0
      ELSE
         IFBGET = -1
      ENDIF
 
C  EXIT ONE WAY OR ANOTHER
C  -----------------------
 
      RETURN
900   CALL ABT('IFBGET - FILE NOT OPEN FOR INPUT')
901   CALL ABT('IFBGET - NO MESSAGE OPEN        ')
      END
