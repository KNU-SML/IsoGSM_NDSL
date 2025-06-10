      SUBROUTINE UFBINX(LUNIT,IMSG,ISUB,USR,I1,I2,IRET,STR)             
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  CUNIT                                               
      CHARACTER*8   SUBSET                                              
      DIMENSION     USR(I1,I2)                                          
      LOGICAL       OPENIT                                              
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      WRITE(CUNIT,'(A8,I2)') 'ON UNIT ',LUNIT                           
                                                                        
C  SEE IF WE NEED TO OPEN A FILE                                        
C  -----------------------------                                        
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      OPENIT = IL.EQ.0                                                  
                                                                        
      IF(OPENIT) THEN                                                   
         CALL OPENBF(LUNIT,'IN',LUNIT)                                  
      ELSE                                                              
         CALL REWNBF(LUNIT,0)                                           
      ENDIF                                                             
                                                                        
C  READ SUBSET ISUB FROM MESSAGE IMSG                                 
C  ------------------------------------                                 
                                                                        
      DO I=1,IMSG-1                                                     
      READ(LUNIT,ERR=900,END=901) MSTR                                  
      ENDDO                                                             
                                                                        
      CALL READMG(LUNIT,SUBSET,JDATE,JRET)                              
      IF(JRET.NE.0) GOTO 901                                            
                                                                        
      DO I=1,ISUB-1                                                     
      IF(NSUB(LUN).GT.MSUB(LUN)) GOTO 902                               
      IBIT = MBYT(LUN)*8                                                
      CALL UPB(NBYT,16,MBAY(1,LUN),IBIT)                                
      MBYT(LUN) = MBYT(LUN) + NBYT                                      
      NSUB(LUN) = NSUB(LUN) + 1                                         
      ENDDO                                                             
      CALL READSB(LUNIT,JRET)                                           
      IF(JRET.NE.0) GOTO 902                                            
                                                                        
      CALL UFBINT(LUNIT,USR,I1,I2,IRET,STR)                             
                                                                        
C  LEAVE THE FILE AS IT WAS BEFORE                                      
C  -------------------------------                                      
                                                                        
      IF(OPENIT) THEN                                                   
         CALL CLOSBF(LUNIT)                                             
      ELSE                                                              
         CALL REWNBF(LUNIT,1)                                           
      ENDIF                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
900   PRINT *,'UFBINX - ERROR READING MESSAGES ', CUNIT
      CALL ABORT
901   PRINT *,'UFBINX - EOF - READING MESSAGES ', CUNIT
      CALL ABORT
902   PRINT *,'UFBINX - EOM - READING SUBSETS  ', CUNIT
      CALL ABORT
      END                                                               
