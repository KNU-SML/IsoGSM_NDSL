      FUNCTION LMSG(SEC0)

      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)

      CHARACTER*8 SEC0

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IMSG = 8/NBYTW
      LMSG = IUPM(SEC0(5:7),24)/8
      IF(LMSG.EQ.0) RETURN
      LMSG = (LMSG+1)*IMSG
      RETURN
      END
