      SUBROUTINE OPENBF(LUNIT,IO,LUNDX)                                 
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TSRIF/  IFIRST                                            
                                                                        
      CHARACTER*(*) IO                                                  
      CHARACTER*20  MSTR                                                
      LOGICAL       SKIPDX,APPEND                                       
                                                                        
      DATA IFIRST/0/                                                    
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(IFIRST.EQ.0) THEN                                              
         CALL WRDLEN                                                    
         CALL BFRINI                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
                                                                        
C  SEE IF A FILE CAN BE OPENED                                          
C  ---------------------------                                          
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(LUN.EQ.0) GOTO 900                                             
      IF(IL .NE.0) GOTO 901                                             
                                                                        
C  CHECK FOR NO BUFR DATA OR NO DATA AT ALL IN AN "IN" FILE             
C  --------------------------------------------------------             
                                                                        
      IF(IO.EQ.'IN' .AND. LUNIT.EQ.LUNDX) THEN                          
         REWIND LUNIT                                                   
         READ(LUNIT,END=100,ERR=902) MSTR                               
         IF(MSTR(1:4).NE.'BUFR') GOTO 902                               
      ENDIF                                                             
                                                                        
C  SET INITIAL OPEN DEFAULTS                                            
C  -------------------------                                            
                                                                        
      REWIND LUNIT                                                      
      NMSG (LUN) = 0                                                    
      NSUB (LUN) = 0                                                    
      MSUB (LUN) = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
      SKIPDX = .FALSE.                                                  
      APPEND = .FALSE.                                                  
                                                                        
C  DECIDE HOW TO SETUP THE DICTIONARY                                   
C  ----------------------------------                                   
                                                                        
      IF(IO.EQ.'IN') THEN                                               
         CALL WTSTAT(LUNIT,LUN,-1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         SKIPDX = LUNIT.NE.LUNDX                                        
      ELSE IF(IO.EQ.'OUT') THEN                                         
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL WRITDX(LUNIT,LUN,LUNDX)                                   
      ELSE IF(IO.EQ.'APN') THEN                                         
         CALL WTSTAT(LUNIT,LUN, 1,0)                                    
         CALL READDX(LUNIT,LUN,LUNDX)                                   
         SKIPDX = .TRUE.                                                
         APPEND = .TRUE.                                                
      ELSE                                                              
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
C  HERE FOR SKIPPING OVER INBEDDED DICTIONARY                           
C  ------------------------------------------                           
                                                                        
      IF(SKIPDX) THEN                                                   
1        READ(LUNIT,END=2) MSTR                                         
         IF(ICHAR(MSTR(17:17)).NE.11) GOTO 2                            
         GOTO 1                                                         
2        BACKSPACE LUNIT                                                
      ENDIF                                                             
                                                                        
C  HERE FOR APPEND MODE                                                 
C  --------------------                                                 
                                                                        
      IF(APPEND) THEN                                                   
3        READ(LUNIT,END=4)                                              
         NMSG(LUN) = NMSG(LUN)+1                                        
         GOTO 3                                                         
4        BACKSPACE LUNIT                                                
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
C  FILE OPENED FOR INPUT IS EMPTY - LET READMG GIVE THE BAD NEWS        
C  -------------------------------------------------------------        
                                                                        
100   REWIND LUNIT                                                      
      CALL WTSTAT(LUNIT,LUN,-1,0)                                       
      CALL DXINIT(LUN,0)                                                
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'OPENBF - TOO MANY FILES OPENED ALREADY       '
      CALL ABORT
901   PRINT *,'OPENBF - FILE ALREADY OPEN                   '
      CALL ABORT
902   PRINT *,'OPENBF - INPUT FILE HAS NON-BUFR DATA        '
      CALL ABORT
903   PRINT *,'OPENBF - IO MUST BE ONE OF "IN" "OUT" "APN"  '
      CALL ABORT
      END                                                               
