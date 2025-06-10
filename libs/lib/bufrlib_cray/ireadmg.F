      FUNCTION IREADMG(LUNIT,SUBSET,IDATE)                              
      CHARACTER*8 SUBSET                                                
      CALL READMG(LUNIT,SUBSET,IDATE,IREADMG)                           
      RETURN                                                            
                                                                        
      ENTRY IREADSB(LUNIT)                                              
      CALL READSB(LUNIT,IREADSB)                                        
      RETURN                                                            
                                                                        
      ENTRY ICOPYSB(LUNIN,LUNOT)                                        
      CALL COPYSB(LUNIN,LUNOT,ICOPYSB)                                  
      RETURN                                                            
                                                                        
      END                                                               
