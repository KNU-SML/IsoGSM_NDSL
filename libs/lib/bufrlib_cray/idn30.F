      FUNCTION IDN30(ADN30,L30)                                         
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*(*) ADN30                                               
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(LEN(ADN30).LT.L30) GOTO 900                                    
      IF(L30.EQ.5) THEN                                                 
         READ(ADN30,'(I5)') IDN30                                       
         IF(IDN30.LT.0 .OR. IDN30.GT.65535) GOTO 901                    
      ELSEIF(L30.EQ.6) THEN                                             
         IDN30 = IFXY(ADN30)                                            
      ELSE                                                              
         GOTO 902                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'IDN30 - FUNCTION INPUT STRING TOO SHORT    '
      CALL ABORT
901   PRINT *,'IDN30 - IDN OUT OF RANGE, NOT A DESCRIPTOR '
      CALL ABORT
902   PRINT *,'IDN30 - CHARACTER LENGTH L30 <> 5 OR 6     '
      CALL ABORT
      END                                                               
