      SUBROUTINE WRDLEN                                                 
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /QUIET / IPRT                                              
                                                                        
      CHARACTER*8 CINT,DINT                                             
      EQUIVALENCE (CINT,INT)                                            
      EQUIVALENCE (DINT,JNT)                                            
      LOGICAL     PRINT                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      PRINT = NBYTW.EQ.0 .AND. IPRT.EQ.1                                
                                                                        
C  COUNT THE BITS IN A WORD - MAX 64 ALLOWED                            
C  -----------------------------------------                            
                                                                        
      INT = 1                                                           
      DO I=1,65                                                         
      INT = ISHFT(INT,1)                                                
      IF(INT.EQ.0) GOTO 10                                              
      ENDDO                                                             
10    IF(I.GE.65)       GOTO 900                                        
      IF(MOD(I,8).NE.0) GOTO 901                                        
      NBITW = I                                                         
      NBYTW = I/8                                                       
                                                                        
C  INDEX THE BYTE STORAGE ORDER -  HIGH BYTE TO LOW BYTE                
C  -----------------------------------------------------                
                                                                        
      JNT = 0                                                           
      DO I=1,NBYTW                                                      
      INT = ISHFT(1,(NBYTW-I)*8)                                        
      DO J=1,NBYTW                                                      
      IF(CINT(J:J).NE.DINT(J:J)) GOTO 20                                
      ENDDO                                                             
20    IF(J.GT.NBYTW) GOTO 902                                           
      IORD(I) = J                                                       
      ENDDO                                                             
                                                                        
C  SET THE NOREVERSE FLAG - 0=NOREVERSE;1=REVERSE                       
C  ----------------------------------------------                       
                                                                        
      NREV = 0                                                          
      DO I=1,NBYTW                                                      
      IF(IORD(I).NE.I) NREV = 1                                         
      ENDDO                                                             
                                                                        
C  SHOW SOME RESULTS                                                    
C  -----------------                                                    
                                                                        
      IF(PRINT) THEN                                                    
         PRINT100,NBYTW,NBITW,NREV,(IORD(I),I=1,NBYTW)                  
      ENDIF                                                             
100   FORMAT(' WRDLEN:NBYTW=',I1,' NBITW=',I2,' IREV=',I1,' IORD=',8I1) 
                                                                        
      RETURN                                                            
900   PRINT *,'WRDLEN - A WORD IS MORE THAN 64 BITS'
      CALL ABORT
901   PRINT *,'WRDLEN - A WORD IS NOT MADE OF BYTES'
      CALL ABORT
902   PRINT *,'WRDLEN - BYTE ORDER CHECKING MISTAKE'
      CALL ABORT
      END                                                               
