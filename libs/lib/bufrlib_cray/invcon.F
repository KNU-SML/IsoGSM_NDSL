      FUNCTION INVCON(NC,LUN,INV1,INV2)                                 
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
                                                                        
      COMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
      TASKCOMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10) 
                                                                        
      REAL*8 VAL,VALS                                                   
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE INVENTORY INTERVAL                                         
C  ----------------------------                                         
                                                                        
      IF(INV1.LE.0 .OR. INV1.GT.NVAL(LUN)) GOTO 99                      
      IF(INV2.LE.0 .OR. INV2.GT.NVAL(LUN)) GOTO 99                      
                                                                        
C  FIND AN OCCURANCE OF NODE IN THE WINDOW MEETING THIS CONDITION       
C  --------------------------------------------------------------       
                                                                        
      DO INVCON=INV1,INV2                                               
      IF(INV(INVCON,LUN).EQ.NODC(NC)) THEN                              
         IF(KONS(NC).EQ.1 .AND. VAL(INVCON,LUN).EQ.VALS(NC)) RETURN     
         IF(KONS(NC).EQ.2 .AND. VAL(INVCON,LUN).NE.VALS(NC)) RETURN     
         IF(KONS(NC).EQ.3 .AND. VAL(INVCON,LUN).LT.VALS(NC)) RETURN     
         IF(KONS(NC).EQ.4 .AND. VAL(INVCON,LUN).GT.VALS(NC)) RETURN     
      ENDIF                                                             
      ENDDO                                                             
                                                                        
99    INVCON = 0                                                        
      RETURN                                                            
      END                                                               
