      SUBROUTINE UPB(NVAL,NBITS,IBAY,IBIT)
 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
 
      DIMENSION IBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      NWD = (IBIT)/NBITW+1
      NBT = MOD(IBIT,NBITW)
      INT = ISHFT(IREV(IBAY(NWD)),NBT)
      INT = ISHFT(INT,NBITS-NBITW)
      LBT = NBT+NBITS
      IF(LBT.GT.NBITW) JNT = IREV(IBAY(NWD+1))
      IF(LBT.GT.NBITW) INT = IOR(INT,ISHFT(JNT,LBT-2*NBITW))
      IBIT = IBIT+NBITS
      NVAL = INT
      RETURN
      END
