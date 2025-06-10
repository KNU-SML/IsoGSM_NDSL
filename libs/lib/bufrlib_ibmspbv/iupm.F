      FUNCTION IUPM(CBAY,NBITS)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      CHARACTER*8 CBAY,CINT
      EQUIVALENCE(CINT,INT)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      IF(NBITS.GT.NBITW) CALL ABT('IUPM - NBITS>WRD LEN')
      CINT = CBAY
      INT  = IREV(INT)
      IUPM = ISHFT(INT,NBITS-NBITW)
      RETURN
      END
