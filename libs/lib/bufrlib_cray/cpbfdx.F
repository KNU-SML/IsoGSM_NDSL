      SUBROUTINE CPBFDX(LUD,LUN)                                        
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB                                                
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION                                    
C  ---------------------------------                                    
                                                                        
      CALL DXINIT(LUN,0)                                                
                                                                        
C  COPY ONE TABLE PARTITION TO ANOTHER                                  
C  -----------------------------------                                  
                                                                        
      INODE(LUN) = INODE(LUD)                                           
                                                                        
      NTBA(LUN) = NTBA(LUD)                                             
      NTBB(LUN) = NTBB(LUD)                                             
      NTBD(LUN) = NTBD(LUD)                                             
                                                                        
      DO I=1,NTBA(LUD)                                                  
      IDNA(I,LUN) = IDNA(I,LUD)                                         
      TABA(I,LUN) = TABA(I,LUD)                                         
      MTAB(I,LUN) = MTAB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBB(LUD)                                                  
      IDNB(I,LUN) = IDNB(I,LUD)                                         
      TABB(I,LUN) = TABB(I,LUD)                                         
      ENDDO                                                             
                                                                        
      DO I=1,NTBD(LUD)                                                  
      IDND(I,LUN) = IDND(I,LUD)                                         
      TABD(I,LUN) = TABD(I,LUD)                                         
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
