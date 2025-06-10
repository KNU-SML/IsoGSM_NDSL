      SUBROUTINE NEMTAB(LUN,NEMO,IDN,TAB,IRET)                          
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*(*) NEMO                                                
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*8   NEMT                                                
      CHARACTER*1   TAB                                                 
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IRET = 0                                                          
      TAB = ' '                                                         
                                                                        
C  LOOK FOR NEMO IN TABLE B                                             
C  ------------------------                                             
                                                                        
      DO I=1,NTBB(LUN)                                                  
      NEMT = TABB(I,LUN)(7:14)                                          
      IF(NEMT.EQ.NEMO) THEN                                             
         IDN  = IDNB(I,LUN)                                             
         TAB  = 'B'                                                     
         IRET = I                                                       
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  LOOK IN TABLE D                                                      
C  ---------------                                                      
                                                                        
      DO I=1,NTBD(LUN)                                                  
      NEMT = TABD(I,LUN)(7:14)                                          
      IF(NEMT.EQ.NEMO) THEN                                             
         IDN  = IDND(I,LUN)                                             
         TAB  = 'D'                                                     
         IRET = I                                                       
         RETURN                                                         
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
