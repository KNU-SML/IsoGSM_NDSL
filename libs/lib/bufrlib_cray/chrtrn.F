      SUBROUTINE CHRTRN(STR,CHR,N)                                      
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*1   CHR(N)                                              
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      DO I=1,N                                                          
      STR(I:I) = CHR(I)                                                 
      ENDDO                                                             
      RETURN                                                            
      END                                                               
