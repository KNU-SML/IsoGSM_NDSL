      SUBROUTINE DATEBF(LUNIT,IY,IM,ID,IH,IDATE)                        
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /TSRIF/  IFIRST                                            
                                                                        
      CHARACTER*24  MSTR                                                
      CHARACTER*8   SEC0                                                
      DIMENSION     MBAY(5000)                                          
      EQUIVALENCE   (SEC0,MBAY(1))                                      
      EQUIVALENCE   (MSTR,MBAY(1))                                      
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
      IF(IFIRST.EQ.0) THEN                                              
         CALL WRDLEN                                                    
         CALL BFRINI                                                    
         IFIRST = 1                                                     
      ENDIF                                                             
                                                                        
      IDATE = -1                                                        
                                                                        
C  SEE IF THE FILE IS ALREADY OPEN TO BUFR INTERFACE (A NO-NO)          
C  -----------------------------------------------------------          
                                                                        
      CALL STATUS(LUNIT,LUN,JL,JM)                                      
      IF(JL.NE.0) THEN
        PRINT *,'DATEBF - FILE ALREADY OPEN'              
        CALL ABORT
      ENDIF
                                                                        
C  READ TO A DATA MESSAGE AND PICK OUT THE DATE                         
C  --------------------------------------------                         
                                                                        
      REWIND LUNIT                                                      
      IMSG = 8/NBYTW+1                                                  
                                                                        
1     READ(LUNIT,ERR=100,END=100) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))      
      IF(MSTR(1:4).NE.'BUFR') GOTO 100                                  
      IF(ICHAR(MSTR(17:17)).EQ.11) GOTO 1                               
                                                                        
C     IY = ICHAR(MSTR(21:21))                                           
      IY = MOD(ICHAR(MSTR(21:21)),100)
      IM = ICHAR(MSTR(22:22))                                           
      ID = ICHAR(MSTR(23:23))                                           
      IH = ICHAR(MSTR(24:24))                                           
      IDATE = IY*1E6+IM*1E4+ID*1E2+IH                                   
                                                                        
      RETURN                                                            
100   CONTINUE
      PRINT *,'DATEBF - THIS IS NOT A BUFR FILE'                    
      CALL ABORT
      END                                                               
