      SUBROUTINE USRTPL(LUN,INVN,NBMP)                                  
                                                                        
      PARAMETER (MAXTMP=5000)                                           
      PARAMETER (MAXINV=90000)                                          
                                                                        
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
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
      DIMENSION    ITMP(MAXTMP),VTMP(MAXTMP)                            
      LOGICAL      DRP,DRS,DRB,DRX                                      
      REAL*8       VAL,VTMP                                             
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     PRINT*,'USRTPL:',LUN,':',INVN,':',NBMP,':',tag(inode(lun))        
                                                                        
      IF(NBMP.LE.0) RETURN                                              
                                                                        
      DRP = .FALSE.                                                     
      DRS = .FALSE.                                                     
      DRX = .FALSE.                                                     
                                                                        
C  SET UP A NODE EXPANSION                                              
C  -----------------------                                              
                                                                        
      IF(INVN.EQ.1) THEN                                                
         NODI = INODE(LUN)                                              
         INV(1,LUN) = NODI                                              
         NVAL(LUN)  = 1                                                 
         IF(NBMP.NE.1) GOTO 900                                         
      ELSEIF(INVN.GT.0 .AND. INVN.LE.NVAL(LUN)) THEN                    
         NODI = INV(INVN,LUN)                                           
         DRP  = TYP(NODI) .EQ. 'DRP'                                    
         DRS  = TYP(NODI) .EQ. 'DRS'                                    
         DRB  = TYP(NODI) .EQ. 'DRB'                                    
         DRX  = DRP .OR. DRS .OR. DRB                                   
         IVAL = VAL(INVN,LUN)                                           
         JVAL = 2**IBT(NODI)-1                                          
         VAL(INVN,LUN) = IVAL+NBMP                                      
         IF(DRB.AND.NBMP.NE.1) GOTO 900                                 
         IF(.NOT.DRX         ) GOTO 901                                 
         IF(IVAL.LT.0.       ) GOTO 902                                 
         IF(IVAL+NBMP.GT.JVAL) GOTO 903                                 
      ELSE                                                              
         GOTO 904                                                       
      ENDIF                                                             
                                                                        
C  RECALL A PRE-FAB NODE EXPANSION SEGMENT                              
C  ---------------------------------------                              
                                                                        
      NEWN = 0                                                          
      N1 = ISEQ(NODI,1)                                                 
      N2 = ISEQ(NODI,2)                                                 
                                                                        
      IF(N1.EQ.0          ) GOTO 905                                    
      IF(N2-N1+1.GT.MAXTMP) GOTO 906                                    
                                                                        
      DO N=N1,N2                                                        
      NEWN = NEWN+1                                                     
      ITMP(NEWN) = JSEQ(N)                                              
      VTMP(NEWN) = VALI(JSEQ(N))                                        
      ENDDO                                                             
                                                                        
C  MOVE OLD NODES - STORE NEW ONES                                      
C  -------------------------------                                      
                                                                        
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) print*,'@:',nval(lun)+newn*nbmp 
      IF(NVAL(LUN)+NEWN*NBMP.GT.MAXINV) GOTO 907                        
                                                                        
CDIR$ IVDEP                                                             
      DO J=NVAL(LUN),INVN+1,-1                                          
      INV(J+NEWN*NBMP,LUN) = INV(J,LUN)                                 
      VAL(J+NEWN*NBMP,LUN) = VAL(J,LUN)                                 
      ENDDO                                                             
                                                                        
      IF(DRP.OR.DRS) VTMP(1) = NEWN                                     
      KNVN = INVN                                                       
                                                                        
      DO I=1,NBMP                                                       
      DO J=1,NEWN                                                       
      KNVN = KNVN+1                                                     
      INV(KNVN,LUN) = ITMP(J)                                           
      VAL(KNVN,LUN) = VTMP(J)                                           
      ENDDO                                                             
      ENDDO                                                             
                                                                        
C  RESET POINTERS AND COUNTERS                                          
C  ---------------------------                                          
                                                                        
      NVAL(LUN) = NVAL(LUN) + NEWN*NBMP                                 
                                                                        
C     print*,tag(inv(invn,lun)),' ',newn,' ',nbmp,' ',nval(lun)         
C     DO I=1,NEWN                                                       
C     PRINT*,TAG(ITMP(I))                                               
C     ENDDO                                                             
                                                                        
                                                                        
      IF(DRX) THEN                                                      
         NODE = NODI                                                    
         INVR = INVN                                                    
4        NODE = JMPB(NODE)                                              
         IF(NODE.GT.0) THEN                                             
            IF(ITP(NODE).EQ.0) THEN                                     
               DO INVR=INVR-1,1,-1                                      
               IF(INV(INVR,LUN).EQ.NODE) THEN                           
                  VAL(INVR,LUN) = VAL(INVR,LUN)+NEWN*NBMP               
                  GOTO 4                                                
               ENDIF                                                    
               ENDDO                                                    
               GOTO 909                                                 
            ELSE                                                        
               GOTO 4                                                   
            ENDIF                                                       
         ENDIF                                                          
      ENDIF                                                             
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'USRTPL - NBMP <> 1 FOR        : ',TAG(NODI)         
      CALL ABORT
901   PRINT *,'USRTPL - NODE NOT SUB,DRP,DRS : ',TAG(NODI)         
      CALL ABORT
902   PRINT *,'USRTPL - NEGATIVE REP FACTOR  : ',TAG(NODI)         
      CALL ABORT
903   PRINT *,'USRTPL - REP FACTOR OVERFLOW  : ',TAG(NODI)         
      CALL ABORT
904   PRINT *,'USRTPL - INVENTORY INDEX OUT OF BOUNDS     '         
      CALL ABORT
905   PRINT *,'USRTPL - UNSET EXPANSION SEG  : ',TAG(NODI)         
      CALL ABORT
906   PRINT *,'USRTPL - TEMP ARRAY OVERFLOW  : ',TAG(NODI)         
      CALL ABORT
907   PRINT *,'USRTPL - INVENTORY OVERFLOW   : ',TAG(NODI)         
      CALL ABORT
908   PRINT *,'USRTPL - TPL CACHE OVERFLOW   : ',TAG(NODI)         
      CALL ABORT
909   PRINT *,'USRTPL - BAD BACKUP STRATEGY  : ',TAG(NODI)         
      CALL ABORT
      END                                                               
