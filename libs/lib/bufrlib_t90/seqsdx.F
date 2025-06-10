      SUBROUTINE SEQSDX(CARD,LUN)                                       
                                                                        
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
                                                                        
      CHARACTER*80  CARD,SEQS                                           
      CHARACTER*12  ATAG,TAGS(50)                                       
      CHARACTER*8   NEMO                                                
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS,TAB                                            
                                                                        
      DATA MAXTGS /50/                                                  
      DATA MAXTAG /12/                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  FIND THE SEQUENCE TAG IN TABLE D AND PARSE THE SEQUENCE STRING       
C  --------------------------------------------------------------       
                                                                        
      NEMO = CARD( 3:10)                                                
      SEQS = CARD(14:78)                                                
                                                                        
      CALL NEMTAB(LUN,NEMO,IDN,TAB,ISEQ)                                
      CALL PARSEQ(SEQS,TAGS,MAXTGS,NTAG)                                
      IF(TAB.NE.'D') GOTO 900                                           
      IF(NTAG.EQ.0 ) GOTO 900                                           
                                                                        
      DO N=1,NTAG                                                       
                                                                        
C  CHECK FOR REPLICATOR AND VALID TAG                                   
C  ----------------------------------                                   
                                                                        
      ATAG = TAGS(N)                                                    
      IREP = 0                                                          
                                                                        
      DO I=1,5                                                          
      IF(ATAG(1:1).EQ.REPS(I,1)) THEN                                   
         DO J=2,MAXTAG                                                  
         IF(ATAG(J:J).EQ.REPS(I,2)) THEN                                
            IF(J.EQ.MAXTAG) GOTO 901                                    
            CALL STRNUM(ATAG(J+1:MAXTAG),NUMR)                          
            IF(I.EQ.1 .AND. NUMR.LE.0  ) GOTO 901                       
            IF(I.EQ.1 .AND. NUMR.GT.255) GOTO 901                       
            IF(I.NE.1 .AND. NUMR.NE.0  ) GOTO 901                       
            ATAG = ATAG(2:J-1)                                          
            IREP = I                                                    
            GOTO 1                                                      
         ENDIF                                                          
         ENDDO                                                          
         GOTO 901                                                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
1     IF(NEMOCK(ATAG).NE.0) GOTO 901                                    
                                                                        
      CALL NEMTAB(LUN,ATAG,IDN,TAB,IRET)                                
      IF(IRET.GT.0) THEN                                                
         IF(TAB.EQ.'B' .AND. IREP.NE.0) GOTO 902                        
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
C  WRITE THE DESCRIPTOR STRING INTO TABD ARRAY                          
C  -------------------------------------------                          
                                                                        
10    IF(IREP.GT.0) CALL PKTDD(ISEQ,LUN,IDNR(IREP,1)+NUMR,IRET)         
      IF(IRET.LT.0) GOTO 904                                            
      CALL PKTDD(ISEQ,LUN,IDN,IRET)                                     
      IF(IRET.LT.0) GOTO 904                                            
                                                                        
      ENDDO                                                             
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'SEQSDX - UNDEFINED SEQUENCE: '             ,   NEMO
      CALL ABORT
901   PRINT *,'SEQSDX - BAD TAG IN SEQUENCE: '            ,TAGS(N)
      CALL ABORT
902   PRINT *,'SEQSDX - REPLICATED ELEMENTS NOT ALLOWED: ',TAGS(N)
      CALL ABORT
903   PRINT *,'SEQSDX - UNDEFINED TAG: '                  ,TAGS(N)
      CALL ABORT
904   PRINT *,'SEQSDX - TOO MANY DESCRIPTORS IN STRING: ' ,NEMO
      CALL ABORT
      END                                                               
