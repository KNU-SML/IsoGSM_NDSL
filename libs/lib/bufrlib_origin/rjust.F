      FUNCTION RJUST(STR)
      CHARACTER*(*) STR
      RJUST = 0
      IF(STR.EQ.' ') RETURN
      LSTR = LEN(STR)
      DO WHILE(STR(LSTR:LSTR).EQ.' ')
         STR = ' '//STR(1:LSTR-1)
      ENDDO
      RETURN
      END
