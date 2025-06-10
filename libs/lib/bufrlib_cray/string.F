      SUBROUTINE STRING(STR,LUN,I1,IO)                                  
                                                                        
      PARAMETER (JCONS=52)                                              
                                                                        
      TASKCOMMON /USRSTR/ JCON(52)                                      
      TASKCOMMON /STCACH/ MSTR,NSTR,LSTR,LUNS(50),USRS(50),ICON(52,50)  
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*80  USRS,UST                                            
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      UST  = STR                                                        
      IF(LEN(STR).GT.80) GOTO 900                                       

C  SEE IF STRING IS IN THE CACHE                                        
C  -----------------------------                                        
                                                                        
      DO I=1,NSTR                                                       
      IF(LUN.EQ.LUNS(I).AND.STR.EQ.USRS(I)) THEN                        
         DO J=1,JCONS                                                   
         JCON(J) = ICON(J,I)                                            
         ENDDO                                                          
         IF(JCON(1).GT.I1) GOTO 901                                     
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  IF NOT PARSE IT AND PUT IT THERE                                     
C  --------------------------------                                     
                                                                        
C     print*,'str=',str                                                 
      CALL PARUSR(STR,LUN,I1,IO)                                        
      LSTR = MAX(MOD(LSTR+1,MSTR+1),1)                                  
      NSTR = MIN(NSTR+1,MSTR)                                           
      LUNS(LSTR) = LUN                                                  
      USRS(LSTR) = STR                                                  
      DO I=1,JCONS                                                      
      ICON(I,LSTR) = JCON(I)                                            
      ENDDO                                                             
      RETURN                                                            
900   PRINT *,'STRING - USER STRING > 80 CHARS:',UST      
      CALL ABORT
901   PRINT *,'STRING - MUST BE AT LEAST I1 STORE NODES:',UST      
      CALL ABORT
      END                                                               
