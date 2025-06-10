      SUBROUTINE UPC(CHR,NCHR,IBAY,IBIT)                                
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*(*) CHR                                                 
      CHARACTER*8   CVAL                                                
      DIMENSION     IBAY(*)                                             
      EQUIVALENCE   (CVAL,IVAL)                                         
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      LB = IORD(NBYTW)                                                  
      DO I=1,NCHR                                                       
      CALL UPB(IVAL,8,IBAY,IBIT)                                        
      CHR(I:I) = CVAL(LB:LB)                                            
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
