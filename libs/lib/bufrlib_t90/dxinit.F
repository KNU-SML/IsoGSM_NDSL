      SUBROUTINE DXINIT(LUN,IOI)                                        
                                                                        
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4                          
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*8   INIB(6,5),INID(5)                                   
      CHARACTER*6   ADN30                                               
      CHARACTER*3   TYPS                                                
      CHARACTER*1   REPS                                                
                                                                        
      DATA INIB   /'------','BYTCNT  ','BYTES  ','+0','+0','16',        
     .             '------','BITPAD  ','NONE   ','+0','+0','1 ',        
     .             '031000','DRF1BIT ','NUMERIC','+0','+0','1 ',        
     .             '031001','DRF8BIT ','NUMERIC','+0','+0','8 ',        
     .             '031002','DRF16BIT','NUMERIC','+0','+0','16'/        
      DATA NINIB  /5/                                                   
                                                                        
      DATA INID   /'        ',                                          
     .             'DRP16BIT',                                          
     .             'DRP8BIT ',                                          
     .             'DRPSTAK ',                                          
     .             'DRP1BIT '/                                          
      DATA NINID  /5/                                                   
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CLEAR OUT A MESSAGE CONTROL WORD PARTITION                           
C  ------------------------------------------                           
                                                                        
      NMSG(LUN)  = 0                                                    
      NSUB(LUN)  = 0                                                    
      MSUB(LUN)  = 0                                                    
      INODE(LUN) = 0                                                    
      IDATE(LUN) = 0                                                    
                                                                        
C  CLEAR OUT A TABLE PARTITION                                          
C  ---------------------------                                          
                                                                        
      NTBA(LUN) = 0                                                     
      DO I=1,NTBA(0)                                                    
      TABA(I,LUN) = ' '                                                 
      MTAB(I,LUN) = 0                                                   
      ENDDO                                                             
                                                                        
      NTBB(LUN) = 0                                                     
      DO I=1,NTBB(0)                                                    
      TABB(I,LUN) = ' '                                                 
      ENDDO                                                             
                                                                        
      NTBD(LUN) = 0                                                     
      DO I=1,NTBD(0)                                                    
      TABD(I,LUN) = ' '                                                 
      CALL PKTDD(I,LUN,0,IRET)                                          
      ENDDO                                                             
                                                                        
      IF(IOI.EQ.0) RETURN                                               
                                                                        
C  INITIALIZE TABLE WITH APRIORI TABLE B AND D ENTRIES                  
C  ---------------------------------------------------                  
                                                                        
      INIB(1,1) = ADN30(IBCT,6)                                         
      INIB(1,2) = ADN30(IPD4,6)                                         
                                                                        
      DO I=1,NINIB                                                      
      NTBB(LUN) = NTBB(LUN)+1                                           
      IDNB(I,LUN) = IFXY(INIB(1,I))                                     
      TABB(I,LUN)(  1:  6) = INIB(1,I)                                  
      TABB(I,LUN)(  7: 70) = INIB(2,I)                                  
      TABB(I,LUN)( 71: 94) = INIB(3,I)                                  
      TABB(I,LUN)( 95: 98) = INIB(4,I)                                  
      TABB(I,LUN)( 99:109) = INIB(5,I)                                  
      TABB(I,LUN)(110:112) = INIB(6,I)                                  
      ENDDO                                                             
                                                                        
      DO I=2,NINID                                                      
      N = NTBD(LUN)+1                                                   
      IDND(N,LUN) = IDNR(I,1)                                           
      TABD(N,LUN)(1: 6) = ADN30(IDNR(I,1),6)                            
      TABD(N,LUN)(7:70) = INID(I)                                       
      CALL PKTDD(N,LUN,IDNR(1,1),IRET)                                  
      CALL PKTDD(N,LUN,IDNR(I,2),IRET)                                  
      NTBD(LUN) = N                                                     
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
