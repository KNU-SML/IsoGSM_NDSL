      SUBROUTINE UPTDD(ID,LUN,IENT,IRET)                                
                                                                        
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LDD = LDXD(IDXV+1)+1                                              
                                                                        
C  CHECK IF IENT IS IN BOUNDS                                           
C  --------------------------                                           
                                                                        
      NDSC = IUPM(TABD(ID,LUN)(LDD:LDD),8)                              
                                                                        
      IF(IENT.EQ.0) THEN                                                
         IRET = NDSC                                                    
         RETURN                                                         
      ELSEIF(IENT.LT.0 .OR. IENT.GT.NDSC) THEN                          
         PRINT *,'UPTDD - IENT OUT OF RANGE'                       
         CALL ABORT
      ENDIF                                                             
                                                                        
C  RETURN THE DESCRIPTOR INDICATED BY IENT                              
C  ---------------------------------------                              
                                                                        
      IDSC = LDD+1 + (IENT-1)*2                                         
      IRET = IUPM(TABD(ID,LUN)(IDSC:IDSC),16)                           
                                                                        
      RETURN                                                            
      END                                                               
