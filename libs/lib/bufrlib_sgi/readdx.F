      SUBROUTINE READDX(LUNIT,LUN,LUNDX)

      COMMON /QUIET/ IPRT

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  GET THE BUFR STATUS OF UNIT LUNDX
C  ---------------------------------

      CALL STATUS(LUNDX,LUD,ILDX,IMDX)

C  READ A DX-TABLE FROM THE INDICATED SOURCE
C  -----------------------------------------

      IF (LUNIT.EQ.LUNDX) THEN
         IF(IPRT.GE.1) PRINT100,LUNDX,LUNIT
         REWIND LUNIT
         CALL RDBFDX(LUNIT,LUN)
      ELSEIF(ILDX.NE.0) THEN
         IF(IPRT.GE.1) PRINT101,LUNDX,LUNIT
         CALL CPBFDX(LUD,LUN)
      ELSEIF(ILDX.EQ.0) THEN
         IF(IPRT.GE.1) PRINT102,LUNDX,LUNIT
         REWIND LUNDX
         CALL RDUSDX(LUNDX,LUN)
      ELSE
         CALL ABORT('READDX - SCREWUP')
      ENDIF

100   FORMAT(' READING BUFR DX-TABLES FROM ',I2,' TO ',I2)
101   FORMAT(' COPYING BUFR DX-TABLES FROM ',I2,' TO ',I2)
102   FORMAT(' READING USER DX-TABLES FROM ',I2,' TO ',I2)

      RETURN
      END
