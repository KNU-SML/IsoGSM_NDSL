      SUBROUTINE RCSTPL(LUN)                                            
                                                                        
      PARAMETER (MAXTMP=1000)                                           
      PARAMETER (MAXRCR=100 )                                           
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      TASKCOMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      TASKCOMMON /USRBIT/ NBIT(INV_),MBIT(INV_)                       
      PARAMETER (MAXINV=INV_)                                          
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION ITMP(MAXTMP,MAXRCR),VTMP(MAXTMP,MAXRCR)                 
      DIMENSION NBMP(2,MAXRCR),NEWN(2,MAXRCR)                           
      DIMENSION KNC(MAXRCR),KNX(MAXRCR)                                 
      REAL*8    VAL,VTMP                                                
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)                                                     
C-----------------------------------------------------------------------
                                                                        
C  SET THE INITIAL VALUES FOR THE TEMPLATE                              
C  ---------------------------------------                              
                                                                        
      INV(1,LUN) = INODE(LUN)                                           
      VAL(1,LUN) = 0                                                    
      NBMP(1,1) = 1                                                     
      NBMP(2,1) = 1                                                     
      NODI = INODE(LUN)                                                 
      NODE = INODE(LUN)                                                 
      MBMP = 1                                                          
      KNVN = 1                                                          
      KNX  = 0                                                          
      NR   = 0                                                          
                                                                        
C  SET UP THE PARAMETRES FOR A LEVEL OF RECURSION                       
C  ----------------------------------------------                       
                                                                        
10    CONTINUE                                                          
                                                                        
      NR = NR+1                                                         
      NBMP(1,NR) = 1                                                    
      NBMP(2,NR) = MBMP                                                 
                                                                        
      N1 = ISEQ(NODE,1)                                                 
      N2 = ISEQ(NODE,2)                                                 
      IF(N1.EQ.0          ) GOTO 905                                    
      IF(N2-N1+1.GT.MAXTMP) GOTO 906                                    
      NEWN(1,NR) = 1                                                    
      NEWN(2,NR) = N2-N1+1                                              
                                                                        
      DO N=1,NEWN(2,NR)                                                 
      NN = JSEQ(N+N1-1)                                                 
      ITMP(N,NR) = NN                                                   
      VTMP(N,NR) = VALI(NN)                                             
      ENDDO                                                             
                                                                        
C  STORE NODES AT SOME RECURSION LEVEL                                  
C  -----------------------------------                                  
                                                                        
20    DO I=NBMP(1,NR),NBMP(2,NR)                                        
      IF(KNX(NR).EQ.0000) KNX(NR) = KNVN                                
      IF(I.GT.NBMP(1,NR)) NEWN(1,NR) = 1                                
      IF(KNVN+NEWN(2,NR)-NEWN(1,NR).GT.MAXINV) GOTO 907                 
      DO J=NEWN(1,NR),NEWN(2,NR)                                        
      KNVN = KNVN+1                                                     
      NODE = ITMP(J,NR)                                                 
      INV(KNVN,LUN) = NODE                                              
      VAL(KNVN,LUN) = VTMP(J,NR)                                        
      MBIT(KNVN) = MBIT(KNVN-1)+NBIT(KNVN-1)                            
      NBIT(KNVN) = IBT(NODE)                                            
      IF(ITP(NODE).EQ.1) THEN                                           
         CALL UPBB(MBMP,NBIT(KNVN),MBIT(KNVN),LUN)                      
         NEWN(1,NR) = J+1                                               
         NBMP(1,NR) = I                                                 
         GOTO 10                                                        
      ENDIF                                                             
      ENDDO                                                             
      NEW = KNVN-KNX(NR)                                                
      VAL(KNX(NR)+1,LUN) = VAL(KNX(NR)+1,LUN) + NEW                     
      KNX(NR) = 0                                                       
      ENDDO                                                             
                                                                        
C  CONTINUE AT ONE RECUSION LEVEL BACK                                  
C  -----------------------------------                                  
                                                                        
      IF(NR-1.NE.0) THEN                                                
         NR = NR-1                                                      
         GOTO 20                                                        
      ENDIF                                                             
                                                                        
C  FINALLY STORE THE LENGTH OF THE SUBSET TEMPLATE                      
C  -----------------------------------------------                      
                                                                        
      NVAL(LUN) = KNVN                                                  
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'RCSTPL - NBMP <> 1 FOR        : ', TAG(NODI)
      CALL ABORT
901   PRINT *,'RCSTPL - NODE NOT SUB,DRP,DRS : ', TAG(NODI)         
      CALL ABORT
902   PRINT *,'RCSTPL - NEGATIVE REP FACTOR  : ', TAG(NODI)         
      CALL ABORT
903   PRINT *,'RCSTPL - REP FACTOR OVERFLOW  : ', TAG(NODI)         
      CALL ABORT
904   PRINT *,'RCSTPL - INVENTORY INDEX OUT OF BOUNDS     '
      CALL ABORT
905   PRINT *,'RCSTPL - UNSET EXPANSION SEG  : ', TAG(NODI)         
      CALL ABORT
906   PRINT *,'RCSTPL - TEMP ARRAY OVERFLOW  : ', TAG(NODI)         
      CALL ABORT
907   PRINT *,'RCSTPL - INVENTORY OVERFLOW   : ', TAG(NODI)         
      CALL ABORT
908   PRINT *,'RCSTPL - TPL CACHE OVERFLOW   : ', TAG(NODI)         
      CALL ABORT
909   PRINT *,'RCSTPL - BAD BACKUP STRATEGY  : ', TAG(NODI)         
      CALL ABORT
      END                                                               
