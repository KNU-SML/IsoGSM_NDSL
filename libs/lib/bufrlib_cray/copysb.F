      SUBROUTINE COPYSB(LUNIN,LUNOT,IRET)                               
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
                                                                        
      CHARACTER*8  SUBSET                                               
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(STATUS,MVB)                                                
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(LUNIN,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      CALL STATUS(LUNOT,LOT,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IL.LT.0) GOTO 904                                              
      IF(IM.EQ.0) GOTO 905                                              
                                                                        
C  MAKE SURE BOTH FILES HAVE THE SAME TABLES                            
C  -----------------------------------------                            
                                                                        
      IF(INODE(LIN).NE.INODE(LOT)) GOTO 906                             
                                                                        
C  SEE IF THERE IS ANOTHER SUBSET IN THE MESSAGE                        
C  ---------------------------------------------                        
                                                                        
      IF(NSUB(LIN).EQ.MSUB(LIN)) THEN                                   
         IRET = -1                                                      
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  COPY THE SUBSET TO THE OUTPUT MESSAGE AND RESET THE POINTERS         
C  ------------------------------------------------------------         
                                                                        
      IBIT = (MBYT(LIN))*8                                              
      CALL UPB(NBYT,16,MBAY(1,LIN),IBIT)                                
                                                                        
      DO N=1,NBYT                                                       
      MBY = MBYT(LIN)+N                                                 
      CALL MVB(MBAY(1,LIN),MBY,IBAY,N)                                  
      ENDDO                                                             
                                                                        
      CALL SUBUPD(LUNOT,LOT,NBYT)                                       
                                                                        
      MBYT(LIN) = MBYT(LIN) + NBYT                                      
      NSUB(LIN) = NSUB(LIN) + 1                                         
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'COPYSB - INPUT FILE IS CLOSED                    '
      CALL ABORT
901   PRINT *,'COPYSB - INPUT FILE IS OPEN FOR OUTPUT           '
      CALL ABORT
902   PRINT *,'COPYSB - NO INPUT FILE MESSAGE OPEN              '
      CALL ABORT
903   PRINT *,'COPYSB - OUTPUT FILE IS CLOSED                   '
      CALL ABORT
904   PRINT *,'COPYSB - OUTPUT FILE IS OPEN FOR INPUT           '
      CALL ABORT
905   PRINT *,'COPYSB - NO OUTPUT FILE MESSAGE OPEN             '
      CALL ABORT
906   PRINT *,'COPYSB - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES'
      CALL ABORT
      END                                                               
