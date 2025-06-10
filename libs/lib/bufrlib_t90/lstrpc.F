      FUNCTION LSTRPC(NODE,LUN)                                         
                                                                        
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
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*3  TYP                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.LT.INODE(LUN) .OR. NODE.GT.ISC(INODE(LUN))) GOTO 900      
                                                                        
      NOD = NODE                                                        
                                                                        
C  FIND THIS OR THE PREVIOUS RPC NODE                                   
C  ----------------------------------                                   
                                                                        
10    IF(TYP(NOD).NE.'RPC') THEN                                        
         NOD = JMPB(NOD)                                                
         IF(NOD.NE.0) GOTO 10                                           
      ENDIF                                                             
                                                                        
      LSTRPC = NOD                                                      
                                                                        
      RETURN                                                            
900   PRINT*,INODE(LUN),':',NODE,':',TAG(NODE)                          
      PRINT *,'LSTRPC - TABLE NODE IS OUT OF BOUNDS'                
      CALL ABORT
      END                                                               
