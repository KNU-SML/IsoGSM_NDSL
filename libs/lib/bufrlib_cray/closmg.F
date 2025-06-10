      SUBROUTINE CLOSMG(LUNIT)                                          
                                                                        
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL MSGWRT(LUNIT,MBAY(1,LUN),MBYT(LUN))              
      CALL WTSTAT(LUNIT,LUN,IL,0)                                       
                                                                        
      RETURN                                                            
900   PRINT *,'CLOSMG - FILE IS CLOSED            '
      CALL ABORT
901   PRINT *,'CLOSMG - FILE IS OPEN FOR INPUT    '
      CALL ABORT
      END                                                               
