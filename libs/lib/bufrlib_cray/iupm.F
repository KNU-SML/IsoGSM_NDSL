      FUNCTION IUPM(CBAY,NBITS)                                         
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 CBAY,CINT                                             
      EQUIVALENCE(CINT,INT)                                             
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NBITS.GT.NBITW) THEN
         PRINT *,'IUPM - NBITS>WRD LEN'
      CALL ABORT
      ENDIF
      CINT = CBAY                                                       
      INT  = IREV(INT)                                                  
      IUPM = ISHFT(INT,NBITS-NBITW)                                     
      RETURN                                                            
      END                                                               
