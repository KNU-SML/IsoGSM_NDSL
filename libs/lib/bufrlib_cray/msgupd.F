      SUBROUTINE MSGUPD(LUNIT,LUN)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      TASKCOMMON /MSGPTR/ NBY0,NBY1,NBY2,NBY3,NBY4,NBY5                 
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
                                                                        
      CHARACTER*8 CBAY                                                  
      EQUIVALENCE (CBAY,JBAY)                                           
                                                                        
C-----------------------------------------------------------------------
      MAXBYT = 9970                                                     
C-----------------------------------------------------------------------
                                                                        
C  PAD THE SUBSET BUFFER                                                
C  ---------------------                                                
                                                                        
      CALL PADIT(IBAY,IBIT,IBYT,8)                                      
      GOTO 1                                                            
                                                                        
C  SPECIAL ENTRY POINT FOR COPYSB                                       
C  ------------------------------                                       
                                                                        
      ENTRY SUBUPD(LUNIT,LUN,JBYT)                                      
      IBYT = JBYT                                                       
                                                                        
C  SEE IF THE NEW SUBSET FITS                                           
C  --------------------------                                           
                                                                        
1     IF(MBYT(LUN)+IBYT.GT.MAXBYT) THEN                                 
         CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))                       
         CALL MSGINI(LUN)                                               
      ENDIF                                                             
                                                                        
      IF(IBYT.GT.MAXBYT-MBYT(LUN)) GOTO 900                             
                                                                        
C  SET A BYTE COUNT AND TRANSFER THE SUBSET BUFFER INTO THE MESSAGE     
C  ----------------------------------------------------------------     
                                                                        
      LBIT = 0                                                          
      CALL PKB(IBYT,16,IBAY,LBIT)                                       
                                                                        
      LBIT = (MBYT(LUN)-4)*8                                            
      DO I=1,IBYT/NBYTW+1                                               
      CALL PKB(IBAY(I),NBITW,MBAY(1,LUN),LBIT)                          
      ENDDO                                                             
                                                                        
C  UPDATE THE SUBSET AND BYTE COUNTERS                                  
C  --------------------------------------                               
                                                                        
      MBYT(LUN)   = MBYT(LUN)   + IBYT                                  
      NSUB(LUN)   = NSUB(LUN)   + 1                                     
                                                                        
      LBIT = (NBY0+NBY1+NBY2+4)*8                                       
      CALL PKB(NSUB(LUN),16,MBAY(1,LUN),LBIT)                           
                                                                        
      LBYT = NBY0+NBY1+NBY2+NBY3                                        
      NBYT = IUPB(MBAY(1,LUN),LBYT+1,24)                                
      LBIT = LBYT*8                                                     
      CALL PKB(NBYT+IBYT,24,MBAY(1,LUN),LBIT)                           
                                                                        
C  RESET THE USER ARRAYS                                                
C  ---------------------                                                
                                                                        
      CALL USRTPL(LUN,1,1)                                              
                                                                        
      RETURN                                                            
900   PRINT *,'MSGUPD - SUBSET LONGER THAN ANY POSSIBLE MESSAGE'
      CALL ABORT
      END                                                               
