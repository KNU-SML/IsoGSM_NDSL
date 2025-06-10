      FUNCTION VALS(STR)

      CHARACTER*(*) STR
      CHARACTER*99  ASTR
      CHARACTER*8   FMT

      data noinline /0/

C----------------------------------------------------------------------
C----------------------------------------------------------------------

      CALL STRSUC(STR,ASTR,LENS)
      WRITE(FMT,'(''(F'',I2,''.0)'')') LENS
      READ(STR,FMT,ERR=900) VAL
      VALS = VAL
      RETURN
900   VALS = 10E10
      RETURN
      END
