      SUBROUTINE NENUCK(NEMO,NUMB,LUN)                                  
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*8   NEMO                                                
      CHARACTER*6   NUMB                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK TABLE A                                                        
C  -------------                                                        
                                                                        
      ENTRY NENUAA(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBA(LUN)                                                  
      IF(NUMB(2:4).EQ.TABA(N,LUN)(1: 3)) GOTO 900                       
      IF(NEMO     .EQ.TABA(N,LUN)(4:11)) GOTO 900                       
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  CHECK TABLE B AND D                                                  
C  -------------------                                                  
                                                                        
      ENTRY NENUBD(NEMO,NUMB,LUN)                                       
                                                                        
      DO N=1,NTBB(LUN)                                                  
      IF(NUMB.EQ.TABB(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABB(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      DO N=1,NTBD(LUN)                                                  
      IF(NUMB.EQ.TABD(N,LUN)(1: 6)) GOTO 900                            
      IF(NEMO.EQ.TABD(N,LUN)(7:14)) GOTO 900                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXIT                                                           
C  ----------                                                           
                                                                        
900   PRINT *,'NENUCK - DUPLICATE NEM/NUM ',NEMO,' ',NUMB        
      CALL ABORT
      END                                                               
