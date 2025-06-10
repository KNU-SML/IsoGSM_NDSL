      SUBROUTINE UFBQCP(LUNIT,QCP,NEMO)                                 
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*6  FXY,ADN30                                            
      CHARACTER*1  TAB                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
                                                                        
      IDN = IFXY('363000')+IFIX(QCP)                                    
      CALL NUMTAB(LUN,IDN,NEMO,TAB,IRET)                                
                                                                        
      RETURN                                                            
900   PRINT *,'UFBQCP - FILE IS CLOSED                       '
      CALL ABORT
      END                                                               
