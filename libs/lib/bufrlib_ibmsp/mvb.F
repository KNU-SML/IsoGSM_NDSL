      SUBROUTINE MVB(IB1,NB1,IB2,NB2,NBM)

      DIMENSION IB1(*),IB2(*)        

C-----------------------------------------------------------------------
CFPP$ EXPAND(UPB,PKB)
C-----------------------------------------------------------------------

      JB1 = 8*(NB1-1)
      JB2 = 8*(NB2-1)

      DO N=1,NBM
      CALL UPB(NVAL,8,IB1,JB1)      
      CALL PKB(NVAL,8,IB2,JB2)      
      ENDDO

      RETURN
      END
