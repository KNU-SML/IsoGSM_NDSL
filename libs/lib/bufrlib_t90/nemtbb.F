      SUBROUTINE NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                   
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*24  UNIT                                                
      CHARACTER*8   NEMO                                                
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(ITAB.LE.0 .OR. ITAB.GT.NTBB(LUN)) GOTO 900                     
                                                                        
C  PULL OUT TABLE B INFORMATION                                         
C  ----------------------------                                         
                                                                        
      IDN  = IDNB(ITAB,LUN)                                             
      NEMO = TABB(ITAB,LUN)( 7:14)                                      
      UNIT = TABB(ITAB,LUN)(71:94)                                      
      ISCL = VAL$(TABB(ITAB,LUN)( 95: 98))                              
      IREF = VAL$(TABB(ITAB,LUN)( 99:109))                              
      IBIT = VAL$(TABB(ITAB,LUN)(110:112))                              
                                                                        
C  CHECK TABLE B CONTENTS                                               
C  ----------------------                                               
                                                                        
      IF(IDN.LT.IFXY('000000')) GOTO 901                                
      IF(IDN.GT.IFXY('063255')) GOTO 901                                
                                                                        
      IF(ISCL.LT.-999 .OR. ISCL.GT.999) GOTO 902                        
      IF(IREF.LE.-1E9 .OR. IREF.GE.1E9) GOTO 903                        
      IF(IBIT.LE.   0 .OR. IBIT.GT. 64) GOTO 904                        
                                                                        
      IF(UNIT.EQ.'CCITT IA5' .AND. MOD(IBIT,8).NE.0) GOTO 905           
                                                                        
      RETURN                                                            
900   PRINT *,'NEMTBB - ITAB NOT IN TABLE B'
      CALL ABORT
901   PRINT *,'NEMTBB - BAD DESCRIPTOR VALUE: ',NEMO               
      CALL ABORT
902   PRINT *,'NEMTBB - BAD SCALE VALUE     : ',NEMO               
      CALL ABORT
903   PRINT *,'NEMTBB - BAD REFERENCE VALUE : ',NEMO               
      CALL ABORT
904   PRINT *,'NEMTBB - BAD BIT WIDTH       : ',NEMO               
      CALL ABORT
905   PRINT *,'NEMTBB - BAD CHAR BIT WIDTH  : ',NEMO               
      CALL ABORT
      END                                                               
