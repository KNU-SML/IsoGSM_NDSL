      SUBROUTINE DRSTPL(INOD,LUN,INV1,INV2,INVN)                        
                                                                        
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
      REAL*8       VAL                                                  
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (INVWIN,USRTPL,NEWWIN)                                     
C-----------------------------------------------------------------------
                                                                        
1     NODE = INOD                                                       
2     NODE = JMPB(NODE)                                                 
C     print*,'drstpl:',tag(inod),':',tag(node),inv1,inv2                
C     print'(5a10)',(tag(inv(i,lun)),i=inv1,inv2)                       
      IF(NODE.EQ.0) GOTO 900                                            
      IF(TYP(NODE).EQ.'DRS' .OR. TYP(NODE).EQ.'DRB') THEN               
         INVN = INVWIN(NODE,LUN,INV1,INV2)                              
         IF(INVN.GT.0) THEN                                             
            CALL USRTPL(LUN,INVN,1)                                     
            CALL NEWWIN(LUN,INV1,INV2)                                  
            INVN = INVWIN(INOD,LUN,INVN,INV2)                           
            IF(INVN.GT.0) RETURN                                        
            GOTO 1                                                      
         ENDIF                                                          
      ENDIF                                                             
      GOTO 2                                                            
900   PRINT *,'DRSTPL - CANT FIND NODE:',TAG(INOD)                 
      CALL ABORT
      END                                                               
