      SUBROUTINE RDBFDX(LUNIT,LUN)                                      
                                                                        
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(50,10),     
     .                IDNA(50,10),IDNB(150,10),IDND(100,10),            
     .                TABA(50,10),TABB(150,10),TABD(100,10)             
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),  
     .                LD30(10),DXSTR(10)                                
                                                                        
      CHARACTER*176 TABD                                                
      CHARACTER*120 TABB,TABB1,TABB2                                    
      CHARACTER*72  TABA                                                
      CHARACTER*56  DXSTR                                               
      CHARACTER*50  DXCMP                                               
      CHARACTER*8   SEC0,NEMO,FMT30                                     
      CHARACTER*6   NUMB,CIDN                                           
      CHARACTER*1   MOCT(24000)                                         
      DIMENSION     MBAY(5000),LDXBD(10),LDXBE(10)                      
      EQUIVALENCE   (MBAY(1),MOCT(1))                                   
                                                                        
      DATA LDXBD /38,70,8*0/                                            
      DATA LDXBE /42,42,8*0/                                            
                                                                        
C-----------------------------------------------------------------------
      JA(I) = IA+1+LDA*(I-1)                                            
      JB(I) = IB+1+LDB*(I-1)                                            
C-----------------------------------------------------------------------
                                                                        
C  INITIALIZE THE DX-TABLE PARTITION AND SOURCE FILE                    
C  -------------------------------------------------                    
                                                                        
      CALL DXINIT(LUN,0)                                                
      REWIND LUNIT                                                      
      IDX = 0                                                           
                                                                        
C  CLEAR THE BUFFER AND READ A MESSAGE                                  
C  -----------------------------------                                  
                                                                        
1     DO I=1,5000                                                       
      MBAY(I) = 0                                                       
      ENDDO                                                             
                                                                        
      IMSG = 8/NBYTW+1                                                  
      READ(LUNIT,ERR=900,END=2) SEC0,(MBAY(I),I=IMSG,LMSG(SEC0))        
      IDX = IDX+1                                                       
                                                                        
C  GET THE SECTION START OCTETS AND LENGTHS                             
C  ----------------------------------------                             
                                                                        
2     I1 = 8                                                            
      L1 = IUPM(MOCT(I1+1),24)                                          
      I2 = I1+L1                                                        
      L2 = IUPM(MOCT(I2+1),24)*IUPM(MOCT(I1+8),1)                       
      I3 = I2+L2                                                        
      L3 = IUPM(MOCT(I3+1),24)                                          
      I4 = I3+L3                                                        
      L4 = IUPM(MOCT(I4+1),24)                                          
                                                                        
C  SEE IF THIS IS A BUFR DX MESSAGE - CHECK FOR RECOGNISABLE DX VERSION 
C  -------------------------------------------------------------------- 
                                                                        
      IF(IUPM(MOCT(I1+9),8).NE.11) THEN                                 
C        BACKSPACE LUNIT                                                
                                                                        
         rewind lunit                                                   
         do ndx=1,idx-1                                                 
         read(lunit,err=910,end=910) sec0,(mbay(i),i=imsg,lmsg(sec0))   
         enddo                                                          
                                                                        
         CALL MAKESTAB                                                  
         RETURN                                                         
      ENDIF                                                             
                                                                        
      IDXS = IUPM(MOCT(I1+10),8)+1                                      
      IF(LDXA(IDXS).EQ.0) GOTO 902                                      
      IF(LDXB(IDXS).EQ.0) GOTO 902                                      
      IF(LDXD(IDXS).EQ.0) GOTO 902                                      
      L30 = LD30(IDXS)                                                  
                                                                        
      DXCMP = ' '                                                       
      CALL CHRTRN(DXCMP,MOCT(I3+8),NXSTR(IDXS))                         
      IF(DXCMP.NE.DXSTR(IDXS)) GOTO 902                                 
                                                                        
C  SECTION 4 - READ DEFINITIONS FOR TABLES A B AND D                    
C  -------------------------------------------------                    
                                                                        
      LDA  = LDXA (IDXS)                                                
      LDB  = LDXB (IDXS)                                                
      LDD  = LDXD (IDXS)                                                
      LDBD = LDXBD(IDXS)                                                
      LDBE = LDXBE(IDXS)                                                
                                                                        
      IA = I4+5                                                         
      LA = IUPM(MOCT(IA),8)                                             
      IB = JA(LA+1)                                                     
      LB = IUPM(MOCT(IB),8)                                             
      ID = JB(LB+1)                                                     
      LD = IUPM(MOCT(ID),8)                                             
                                                                        
C  TABLE A                                                              
C  -------                                                              
                                                                        
      DO I=1,LA                                                         
      N = NTBA(LUN)+1                                                   
      IF(N.GT.NTBA(0)) GOTO 903                                         
      CALL CHRTRN(TABA(N,LUN),MOCT(JA(I)),LDA)                          
      NUMB = TABA(N,LUN)(1:3)                                           
      NEMO = TABA(N,LUN)(4:11)                                          
      CALL NENUAA(NEMO,NUMB,LUN)                                        
      READ(NUMB(1:3),'(I3)') IDNA(N,LUN)                                
      NTBA(LUN) = N                                                     
      ENDDO                                                             
                                                                        
C  TABLE B                                                              
C  -------                                                              
                                                                        
      DO I=1,LB                                                         
      N = NTBB(LUN)+1                                                   
      IF(N.GT.NTBB(0)) GOTO 904                                         
      CALL CHRTRN(TABB1,MOCT(JB(I)     ),LDBD)                          
      CALL CHRTRN(TABB2,MOCT(JB(I)+LDBD),LDBE)                          
      TABB(N,LUN) = TABB1(1:LDXBD(IDXV+1))//TABB2(1:LDXBE(IDXV+1))      
      NUMB = TABB(N,LUN)(1:6)                                           
      NEMO = TABB(N,LUN)(7:14)                                          
      CALL NENUBD(NEMO,NUMB,LUN)                                        
      IDNB(N,LUN) = IFXY(NUMB)                                          
      NTBB(LUN) = N                                                     
      ENDDO                                                             
                                                                        
C  TABLE D                                                              
C  -------                                                              
                                                                        
      DO I=1,LD                                                         
      N = NTBD(LUN)+1                                                   
      IF(N.GT.NTBD(0)) GOTO 905                                         
      CALL CHRTRN(TABD(N,LUN),MOCT(ID+1),LDD)                           
      NUMB = TABD(N,LUN)(1:6)                                           
      NEMO = TABD(N,LUN)(7:14)                                          
      CALL NENUBD(NEMO,NUMB,LUN)                                        
      IDND(N,LUN) = IFXY(NUMB)                                          
      ND = IUPM(MOCT(ID+LDD+1),8)                                       
      IF(ND.GT.50) GOTO 906                                             
      DO J=1,ND                                                         
      NDD = ID+LDD+2 + (J-1)*L30                                        
      CALL CHRTRN(CIDN,MOCT(NDD),L30)                                   
      IDN = IDN30(CIDN,L30)                                             
      CALL PKTDD(N,LUN,IDN,IRET)                                        
      IF(IRET.LT.0) GOTO 908                                            
      ENDDO                                                             
      ID = ID+LDD+1 + ND*L30                                            
      IF(IUPM(MOCT(ID+1),8).EQ.0) ID = ID+1                             
      NTBD(LUN) = N                                                     
      ENDDO                                                             
                                                                        
C  GOTO READ THE NEXT MESSAGE                                           
C  --------------------------                                           
                                                                        
      GOTO 1                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'RDBFDX - I/O ERROR READING DX MESSAGE          '
      CALL ABORT
901   PRINT *,'RDBFDX - EOF >>>>> READING DX MESSAGE          '
      CALL ABORT
902   PRINT *,'RDBFDX - UNEXPECTED DX MESSAGE TYPE OR CONTENTS'
      CALL ABORT
903   PRINT *,'RDBFDX - TOO MANY TABLE A ENTRIES              '
      CALL ABORT
904   PRINT *,'RDBFDX - TOO MANY TABLE B ENTRIES              '
      CALL ABORT
905   PRINT *,'RDBFDX - TOO MANY TABLE D ENTRIES              '
      CALL ABORT
906   PRINT *,'RDBFDX - TOO MANY DESCRIPTORS IN TABLE D ENTRY '
      CALL ABORT
907   PRINT *,'RDBFDX - ERROR READING IDN SEQ FROM MOCT       '
      CALL ABORT
908   PRINT *,'RDBFDX - BAD RETURN FROM PKTDD                 '
      CALL ABORT
909   PRINT *,'RDBFDX - DESC COUNT IN TABD <> MOCT            '
      CALL ABORT
910   PRINT *,'RDBFDX - ERR/EOF POSITIONING AFTER DX MESSAGES '
      CALL ABORT
      END                                                               
