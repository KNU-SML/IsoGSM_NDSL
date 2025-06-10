      FUNCTION IUPB(MBAY,NBYT,NBIT)
 
      DIMENSION MBAY(*)
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
      MBIT = (NBYT-1)*8
      CALL UPB(IUPB,NBIT,MBAY,MBIT)
      RETURN
      END
