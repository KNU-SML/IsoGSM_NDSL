      SUBROUTINE NEMTBD(LUN,ITAB,NSEQ,NEMS,IRPS,KNTS)                   
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*8   NEMO,NEMS,NEMT                                      
      CHARACTER*1   TAB                                                 
      DIMENSION     NEMS(50),IRPS(50),KNTS(50)                          
      LOGICAL       REP                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(ITAB.LE.0 .OR. ITAB.GT.NTBD(LUN)) GOTO 900                     
                                                                        
      REP  = .FALSE.                                                    
                                                                        
C  CLEAR THE RETURN VALUES                                              
C  -----------------------                                              
                                                                        
      NSEQ = 0                                                          
                                                                        
      DO I=1,50                                                         
      NEMS(I) = ' '                                                     
      IRPS(I) = 0                                                       
      KNTS(I) = 0                                                       
      ENDDO                                                             
                                                                        
C  PARSE THE TABLE D ENTRY                                              
C  -----------------------                                              
                                                                        
      NEMO = TABD(ITAB,LUN)(7:14)                                       
      IDSC = IDND(ITAB,LUN)                                             
      CALL UPTDD(ITAB,LUN,0,NDSC)                                       
                                                                        
      IF(IDSC.LT.IFXY('300000')) GOTO 901                               
      IF(IDSC.GT.IFXY('363255')) GOTO 901                               
C     IF(NDSC.LE.0             ) GOTO 902                               
                                                                        
      DO J=1,NDSC                                                       
      IF(NSEQ+1.GT.50) GOTO 903                                         
      CALL UPTDD(ITAB,LUN,J,IDSC)                                       
      CALL NUMTAB(LUN,IDSC,NEMT,TAB,IRET)                               
      IF(TAB.EQ.'R') THEN                                               
         IF(REP) GOTO 904                                               
         REP = .TRUE.                                                   
         IF(IRET.LT.0) THEN                                             
            IRPS(NSEQ+1) = 1                                            
            KNTS(NSEQ+1) = ABS(IRET)                                    
         ELSEIF(IRET.GT.0) THEN                                         
            IRPS(NSEQ+1) = IRET                                         
         ENDIF                                                          
      ELSEIF(TAB.EQ.'F') THEN                                           
         IF(.NOT.REP) GOTO 904                                          
         IRPS(NSEQ+1) = IRET                                            
         REP = .FALSE.                                                  
      ELSEIF(TAB.EQ.'B' .OR. TAB.EQ.'D') THEN                           
         REP = .FALSE.                                                  
         NSEQ = NSEQ+1                                                  
         NEMS(NSEQ) = NEMT                                              
      ELSE                                                              
         GOTO 905                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'NEMTBD - ITAB NOT IN TABLE D   '  
      CALL ABORT
901   PRINT *,'NEMTBD - BAD DESCRIPTOR VALUE: ' ,NEMO    
      CALL ABORT
902   PRINT *,'NEMTBD - ZERO LENGTH SEQUENCE: '          ,NEMO     
      CALL ABORT
903   PRINT *,'NEMTBD - TOO MANY DESCRIPTORS IN SEQ: '   ,NEMO     
      CALL ABORT
904   PRINT *,'NEMTBD - REPLICATOR OUT OF ORDER IN SEQ: ',NEMO     
      CALL ABORT
905   PRINT *,'NEMTBD - BAD DESCRIPTOR IN SEQUENCE: '    ,NEMO     
      CALL ABORT
      END                                                              
