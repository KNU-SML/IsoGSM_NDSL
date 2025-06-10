      SUBROUTINE WRITSB(LUNIT)                                          
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.EQ.0) GOTO 902                                              
                                                                        
C  PACK UP THE SUBSET AND PUT IT INTO THE MESSAGE                       
C  ----------------------------------------------                       
                                                                        
      CALL WRTREE(LUN)                                                  
      CALL MSGUPD(LUNIT,LUN)                                            
                                                                        
      RETURN                                                            
900   PRINT *,'WRITSB - FILE IS CLOSED                     '
      CALL ABORT
901   PRINT *,'WRITSB - FILE IS OPEN FOR INPUT             '
      CALL ABORT
902   PRINT *,'WRITSB - NO MESSAGE OPEN                    '
      CALL ABORT
      END                                                               
