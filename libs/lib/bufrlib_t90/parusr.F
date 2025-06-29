      SUBROUTINE PARUSR(STR,LUN,I1,IO)                                  
                                                                        
      TASKCOMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10) 
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*80  UST                                                 
      CHARACTER*20  UTG(30)                                             
      LOGICAL       BUMP                                                
      REAL*8        VALS                                                
                                                                        
      DATA MAXUSR /30/                                                  
      DATA MAXNOD /20/                                                  
      DATA MAXCON /10/                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      UST  = STR
      IF(LEN(STR).GT.80) GOTO 900                                       
                                                                        
      NCON = 0                                                          
      NNOD = 0                                                          
                                                                        
C  PROCESS STRING PIECES(S) INTO COND NODES AND STORE NODES             
C  --------------------------------------------------------             
                                                                        
      CALL PARSEQ(UST,UTG,MAXUSR,NTOT)                                  
                                                                        
      DO N=1,NTOT                                                       
      CALL PARUTG(LUN,IO,UTG(N),NOD,KON,VAL,*908)                       
      IF(KON.NE.0) THEN                                                 
         NCON = NCON+1                                                  
         IF(NCON.GT.MAXCON) GOTO 901                                    
         NODC(NCON) = NOD                                               
         KONS(NCON) = KON                                               
         VALS(NCON) = VAL                                               
      ELSE                                                              
         NNOD = NNOD+1                                                  
         IF(NNOD.GT.MAXNOD) GOTO 902                                    
         NODS(NNOD) = NOD                                               
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  SORT COND NODES IN JUMP/LINK TABLE ORDER                             
C  ----------------------------------------                             
                                                                        
      DO I=1,NCON                                                       
      DO J=I+1,NCON                                                     
      IF(NODC(I).GT.NODC(J)) THEN                                       
         NOD     = NODC(I)                                              
         NODC(I) = NODC(J)                                              
         NODC(J) = NOD                                                  
                                                                        
         KON     = KONS(I)                                              
         KONS(I) = KONS(J)                                              
         KONS(J) = KON                                                  
                                                                        
         VAL     = VALS(I)                                              
         VALS(I) = VALS(J)                                              
         VALS(J) = VAL                                                  
      ENDIF                                                             
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  CHECK ON SPECIAL RULES FOR BUMP NODES                                
C  -------------------------------------                                
                                                                        
      BUMP = .FALSE.                                                    
                                                                        
      DO N=1,NCON                                                       
      IF(KONS(N).EQ.5) THEN                                             
         IF(IO.EQ.0)   GOTO 903                                         
         IF(N.NE.NCON) GOTO 904                                         
         BUMP = .TRUE.                                                  
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  CHECK STORE NODE COUNT AND ALIGNMENT                                 
C  ------------------------------------                                 
                                                                        
      IF(.NOT.BUMP .AND. NNOD.EQ.0) GOTO 905                            
      IF(NNOD.GT.I1)                GOTO 906                            
                                                                        
      IRPC = -1                                                         
      DO I=1,NNOD                                                       
      IF(NODS(I).GT.0) THEN                                             
         IF(IRPC.LT.0) IRPC = LSTRPC(NODS(I),LUN)                       
         IF(IRPC.NE.LSTRPC(NODS(I),LUN)) GOTO 907                       
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'PARUSR - USER STRING > 80 CHARS         :',UST      
      CALL ABORT
901   PRINT *,'PARUSR - TOO MANY COND NODES            :',UST
      CALL ABORT
902   PRINT *,'PARUSR - TOO MANY STOR NODES            :',UST      
      CALL ABORT
903   PRINT *,'PARUSR - BUMP ON INPUT NOT ALLOWED      :',UST      
      CALL ABORT
904   PRINT *,'PARUSR - BUMP MUST BE ON INNER NODE     :',UST      
      CALL ABORT
905   PRINT *,'PARUSR - USER STRING HAS NO STORE NODES :',UST      
      CALL ABORT
906   PRINT *,'PARUSR - MUST BE AT LEAST I1 STORE NODES:',UST      
      CALL ABORT
907   PRINT *,'PARUSR - STORE NODES MUST IN ONE REP GRP:',UST      
      CALL ABORT
908   PRINT *,'PARUSR - PARUTG:'                         ,UST      
      CALL ABORT
      END                                                               
