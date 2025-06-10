      SUBROUTINE WTSTAT(LUNIT,LUN,IL,IM)

      COMMON /STBFR/ IOLUN(10),IOMSG(10)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK ON THE ARGUMENTS
C  ----------------------

      IF(LUNIT.LE.0)            GOTO 900
      IF(LUN  .LE.0)            GOTO 901
      IF(IL.LT.-1 .OR. IL.GT.1) GOTO 902
      IF(IM.LT. 0 .OR. IL.GT.1) GOTO 903

C  CHECK ON LUNIT-LUN COMBINATION
C  ------------------------------

      IF(ABS(IOLUN(LUN)).NE.LUNIT) THEN
         IF(IOLUN(LUN).NE.0) GOTO 905
      ENDIF

C  RESET THE FILE STATUSES
C  -----------------------

      IF(IL.NE.0) THEN
         IOLUN(LUN) = SIGN(LUNIT,IL)
         IOMSG(LUN) = IM
      ELSE
         IOLUN(LUN) = 0
         IOMSG(LUN) = 0
      ENDIF

      RETURN
900   CALL ABT('WTSTAT - BAD LUNIT                               ')
901   CALL ABT('WTSTAT - BAD LUN                                 ')
902   CALL ABT('WTSTAT - BAD IL                                  ')
903   CALL ABT('WTSTAT - BAD IM                                  ')
905   CALL ABT('WTSTAT - ATTEMPT TO REDEFINE EXISITING FILE UNIT ')
      END
