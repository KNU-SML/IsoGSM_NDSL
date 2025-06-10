      SUBROUTINE CHEKSTAB(LUN)                                          
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO,NEMS(50)                                       
      CHARACTER*1   TAB                                                 
      DIMENSION     IRPS(50),KNTS(50)                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  THERE MUST BE ENTRIES IN TABLES A, B, AND D                          
C  -------------------------------------------                          
                                                                        
      IF(NTBA(LUN).EQ.0) GOTO 900                                       
      IF(NTBB(LUN).EQ.0) GOTO 901                                       
      IF(NTBD(LUN).EQ.0) GOTO 902                                       
                                                                        
C  MAKE SURE EACH TABLE A ENTRY DEFINED AS A SEQUENCE                   
C  --------------------------------------------------                   
                                                                        
      DO I=1,NTBA(LUN)                                                  
      NEMO = TABA(I,LUN)(4:11)                                          
      CALL NEMTAB(LUN,NEMO,IDN,TAB,IRET)                                
      IF(TAB.NE.'D') GOTO 903                                           
      ENDDO                                                             
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBB(LUN)                                               
      CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                         
      ENDDO                                                             
                                                                        
C  CHECK TABLE D CONTNETS                                               
C  ----------------------                                               
                                                                        
      DO ITAB=1,NTBD(LUN)                                               
      CALL NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'CHEKSTAB - EMPTY TABLE A'
      CALL ABORT
901   PRINT *,'CHEKSTAB - EMPTY TABLE B'
      CALL ABORT
902   PRINT *,'CHEKSTAB - EMPTY TABLE D'
      CALL ABORT
903   PRINT *,'CHEKSTAB - NO SEQUENCE DEFINED FOR TABLE A: ',NEMO
      CALL ABORT
      END                                                               
