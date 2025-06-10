      SUBROUTINE UFBPKS(LUBFR,LUBFO)                                    
                                                                        
      PARAMETER(MAXMSG=100)                                             
                                                                        
      COMMON /QUIET/IPRT                                                
                                                                        
      CHARACTER*8 SUBSET,MSGTYP(MAXMSG)                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IPRT = 1                                                          
                                                                        
C  OPEN THE INPUT AND OUTPUT FILES                                      
C  -------------------------------                                      
                                                                        
      CALL OPENBF(LUBFR,'IN ',LUBFR)                                    
      CALL OPENBF(LUBFO,'OUT',LUBFR)                                    
                                                                        
C  MAKE THE MESSAGE TYPE TABLE                                          
C  ---------------------------                                          
                                                                        
      NMSG = 0                                                          
                                                                        
10    DO WHILE(IREADMG(LUBFR,SUBSET,IDATE).EQ.0)                        
      DO N=1,NMSG                                                       
      IF(SUBSET.EQ.MSGTYP(N)) GOTO 10                                   
      ENDDO                                                             
      IF(NMSG+1.GT.MAXMSG) GOTO 900                                     
      NMSG = NMSG+1                                                     
      MSGTYP(NMSG) = SUBSET                                             
      ENDDO                                                             
      CALL CLOSBF(LUBFR)                                                
                                                                        
C  FOR EACH MESSAGE TYPE CONSOLIDATE ALL OCCURANCES INTO LUBFO          
C  -----------------------------------------------------------          
                                                                        
      IPRT = 0                                                          
                                                                        
      DO N=1,NMSG                                                       
      PRINT*,'PACKING MSGTYP:',MSGTYP(N)                                
      CALL OPENBF(LUBFR,'IN',LUBFO)                                     
      DO WHILE(IREADMG(LUBFR,SUBSET,IDATE).EQ.0)                        
      IF(SUBSET.EQ.MSGTYP(N)) THEN                                      
         CALL OPENMB(LUBFO,SUBSET,IDATE)                                
         DO WHILE(ICOPYSB(LUBFR,LUBFO).EQ.0)                            
         ENDDO                                                          
      ENDIF                                                             
      ENDDO                                                             
      CALL CLOSBF(LUBFR)                                                
      ENDDO                                                             
      CALL CLOSBF(LUBFO)                                                
                                                                        
C  EXIT(S)                                                              
C  -------                                                              
                                                                        
      RETURN                                                            
900   PRINT *,'UFBPKS - NMSG>MAXMSG'
      CALL ABORT
      END                                                               
