      SUBROUTINE OPENMG(LUNIT,SUBSET,JDATE)                             
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
                                                                        
      CHARACTER*(*) SUBSET                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUS                                                
C  ---------------------                                                
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.LT.0) GOTO 901                                              
      IF(IM.NE.0) CALL CLOSMG(LUNIT)                                    
      CALL WTSTAT(LUNIT,LUN,IL, 1)                                      
                                                                        
C  GET SOME SUBSET PARTICULARS                                          
C  ---------------------------                                          
                                                                        
      CALL NEMTBA(LUN,SUBSET,MTYP,INOD)                                 
      IF(INODE(LUN).NE.INOD) CALL STRCLN                                
      INODE(LUN) = INOD                                                 
      IDATE(LUN) = JDATE                                                
                                                                        
C  INITIALIZE THE OPEN MESSAGE                                          
C  ---------------------------                                          
                                                                        
      CALL MSGINI(LUN)                                                  
      CALL USRTPL(LUN,1,1)                                              
                                                                        
      RETURN                                                            
900   PRINT *,'OPENMG - FILE IS CLOSED            '
      CALL ABORT
901   PRINT *,'OPENMG - FILE IS OPEN FOR INPUT    '
      CALL ABORT
      END                                                               
