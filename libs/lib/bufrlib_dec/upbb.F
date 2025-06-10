      SUBROUTINE UPBB(INT,NBIT,MBIT,LUN)
 
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      NWD = MBIT/NBITW + 1
      NBT = MOD(MBIT,NBITW)
      LBT = NBT+NBIT
      IBA = IREV(MBAY(NWD,LUN))
      INT = ISHFT(ISHFT(IBA,NBT),NBIT-NBITW)
      IF(LBT.GT.NBITW) THEN
         IBA = IREV(MBAY(NWD+1,LUN))
         INT = IOR(INT,ISHFT(IBA,LBT-2*NBITW))
      ENDIF
      RETURN
      END
