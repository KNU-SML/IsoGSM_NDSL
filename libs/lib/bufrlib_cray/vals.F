      FUNCTION VAL$(STR)                                                
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*99  ASTR                                                
      CHARACTER*8   FMT                                                 
                                                                        
      data noinline /0/                                                 
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      ASTR = STR                                                        
      WRITE(FMT,'(''(F'',I2,''.0)'')') LEN(STR)                         
      READ(STR,FMT,ERR=900) VAL                                         
      VAL$ = VAL                                                        
      RETURN                                                            
900   VAL$ = 10E10                                                      
      RETURN                                                            
      END                                                               
