      SUBROUTINE NXTWIN(LUN,IWIN,JWIN)                                  
                                                                        
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      TASKCOMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
cfpp$ expand (lstrpc)                                                   
C---------------------------------------------------------------------- 
                                                                        
      IF(JWIN.EQ.NVAL(LUN)) THEN                                        
         IWIN = 0                                                       
         RETURN                                                         
      ENDIF                                                             
                                                                        
C  FIND THE NEXT SEQUENTIAL WINDOW                                      
C  -------------------------------                                      
                                                                        
      NODE = INV(IWIN,LUN)                                              
      IF(LSTRPC(NODE,LUN).NE.NODE) print*,'bad node=',node,iwin         
      IF(LSTRPC(NODE,LUN).NE.NODE) THEN
      PRINT *,'NXTWIN - NOT RPC'       
      CALL ABORT
      ENDIF
      IF(VAL(JWIN,LUN).EQ.0) THEN                                       
         IWIN = 0                                                       
      ELSE                                                              
         IWIN = JWIN                                                    
         JWIN = IWIN+VAL(IWIN,LUN)                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
