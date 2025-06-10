      SUBROUTINE UFBINT(LUNIO,USR,I1,I2,IRET,STR)                       
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      TASKCOMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
                                                                        
      CHARACTER*(*) STR                                                 
      DIMENSION     USR(I1,I2)                                          
      REAL*8        USR,VAL                                             
                                                                        
C---------------------------------------------------------------------- 
CFPP$ EXPAND (STATUS,UFBRW)                                             
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      LUNIT = ABS(LUNIO)                                                
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
      IO = MIN(MAX(0,IL),1)                                             
      IF(LUNIO.NE.LUNIT) IO = 0                                         
                                                                        
C  PARSE THE INPUT STRING - READ/WRITE VALUES                           
C  ------------------------------------------                           
                                                                        
      CALL STRING(STR,LUN,I1,IO)                                        
      CALL UFBRW(LUN,USR,I1,I2,IO,IRET)                                 
                                                                        
      IF(IO.EQ.1 .AND. IRET.NE.I2) THEN                                 
         PRINT*,STR                                                     
         GOTO 903                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   CONTINUE
      PRINT *,'UFBINT - FILE IS CLOSED'        
      CALL ABORT
901   PRINT *,'UFBINT - NO MESSAGE OPEN'        
      CALL ABORT
902   PRINT *,'UFBINT - I-NODE MISMATCH'
      CALL ABORT
903   PRINT *,'UFBINT - INCOMPLETE WRITE'        
      CALL ABORT
      END                                                               
