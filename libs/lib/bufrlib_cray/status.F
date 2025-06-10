      SUBROUTINE STATUS(LUNIT,LUN,IL,IM)                                
                                                                        
      COMMON /STBFR/ IOLUN(10),IOMSG(10)                                
      TASKCOMMON /STBFM/ MOLUN(10),MOMSG(10)                            
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(LUNIT.LE.0 .OR. LUNIT.GT.99) GOTO 900                          
                                                                        
C  CLEAR THE STATUS INDICATORS                                          
C  ---------------------------                                          
                                                                        
      LUN = 0                                                           
      IL  = 0                                                           
      IM  = 0                                                           
                                                                        
C  SEE IF THE UNIT IS DEFINED                                           
C  --------------------------                                           
                                                                        
      DO I=1,10                                                         
      IF(ABS(IOLUN(I)).EQ.LUNIT) LUN = I                                
      ENDDO                                                             
                                                                        
C  IF NOT, CHECK FOR FILE SPACE - RETURN LUN=0 IF NO FILE SPACE         
C  ------------------------------------------------------------         
                                                                        
      IF(LUN.EQ.0) THEN                                                 
         DO I=1,10                                                      
         IF(IOLUN(I).EQ.0) THEN                                         
            LUN = I                                                     
            RETURN                                                      
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  IF FILE DEFINED RETURN STATUSES                                      
C  -------------------------------                                      
                                                                        
      IL = SIGN(1,IOLUN(LUN))                                           
      IM = MOMSG(LUN)                                                   
                                                                        
      RETURN                                                            
900   PRINT *,'STATUS - ILLEGAL UNIT GIVEN'
      CALL ABORT
      END                                                               
