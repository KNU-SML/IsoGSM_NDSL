      SUBROUTINE DATEBF(LUNIT,IY,IM,ID,IH,IDATE)

      CHARACTER*24  MSTR

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IDATE = -1

C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)
C  -----------------------------------------------------------

      CALL STATUS(LUNIT,LUN,JL,JM)
      IF(JL.NE.0) CALL ABORT('DATEBF - FILE ALREADY OPEN')

C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL
C  ----------------------------------------

      REWIND LUNIT
      READ(LUNIT,END=100,ERR=100) MSTR
      IF(MSTR(1:4).NE.'BUFR') GOTO 100

C  READ TO A DATA MESSAGE AND PICK OUT THE DATE
C  --------------------------------------------

1     READ(LUNIT,END=100,ERR=100) MSTR
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1
      IY = ICHAR(MSTR(21:21))
      IY = MOD(ICHAR(MSTR(21:21)),100)
      IM = ICHAR(MSTR(22:22))
      ID = ICHAR(MSTR(23:23))
      IH = ICHAR(MSTR(24:24))
      IDATE = IY*1E6+IM*1E4+ID*1E2+IH
100   RETURN
      END
