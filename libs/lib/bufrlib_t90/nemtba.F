      SUBROUTINE NEMTBA(LUN,NEMO,MTYP,INOD)                             
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*20  NEMT                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      NEMT = NEMO                                                       
      IRET = 0                                                          
                                                                        
C  LOOK FOR NEMO IN TABLE A                                             
C  ------------------------                                             
                                                                        
      DO I=1,NTBA(LUN)                                                  
      IF(TABA(I,LUN)(4:11).EQ.NEMO) THEN                                
         MTYP = IDNA(I,LUN)                                             
         INOD = MTAB(I,LUN)                                             
         IF(MTYP.LE.0 .OR. MTYP.GT.999) GOTO 900                        
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      PRINT *,'NEMTBA - CANT FIND ',NEMT
      CALL ABORT
900   PRINT *,'NEMTBA - BAD MTYP  ',NEMT
      CALL ABORT
      END                                                               
