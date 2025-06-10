      FUNCTION LMSG(SEC0)                                               
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
                                                                        
      CHARACTER*8 SEC0                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      LMSG = IUPM(SEC0(5:7),24)/NBYTW+1                                 
      RETURN                                                            
      END                                                               
