      SUBROUTINE IPKM(CBAY,NBYT,N)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 CBAY,CINT                                             
      EQUIVALENCE(CINT,INT)                                             
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NBYT.GT.NBYTW) THEN
         PRINT *,'IPKM - NBYT>WRD LEN'
      CALL ABORT
      ENDIF
      INT = IREV(ISHFT(N,(NBYTW-NBYT)*8))                               
      DO I=1,NBYT                                                       
      CBAY(I:I) = CINT(I:I)                                             
      ENDDO                                                             
      RETURN                                                            
      END                                                               
