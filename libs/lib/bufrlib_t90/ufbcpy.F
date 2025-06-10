      SUBROUTINE UFBCPY(LUBIN,LUBOT)                                    
                                                                        
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      TASKCOMMON /USRINT/ NVAL(10),INV(INV_,10),VAL(INV_,10)
                                                                        
      REAL*8 VAL                                                        
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUSES AND I-NODE                                   
C  ----------------------------------                                   
                                                                        
      CALL STATUS(LUBIN,LUI,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUI).NE.INV(1,LUI)) GOTO 902                             
                                                                        
      CALL STATUS(LUBOT,LUO,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IM.EQ.0) GOTO 904                                              
      IF(INODE(LUI).NE.INODE(LUO)) GOTO 905                             
                                                                        
C  EVERYTHING OKAY COPY USER ARRAY FROM LUI TO LUO                      
C  -----------------------------------------------                      
                                                                        
      NVAL(LUO) = NVAL(LUI)                                             
                                                                        
      DO N=1,NVAL(LUI)                                                  
      INV(N,LUO) = INV(N,LUI)                                           
      VAL(N,LUO) = VAL(N,LUI)                                           
      ENDDO                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'UFBCPY - INPUT  FILE IS NOT OPEN'
      CALL ABORT
901   PRINT *,'UFBCPY - INPUT  MESG IS NOT OPEN'
      CALL ABORT
902   PRINT *,'UFBCPY - INPUT  I-NODE  MISMATCH'
      CALL ABORT
903   PRINT *,'UFBCPY - OUTPUT FILE IS NOT OPEN'
      CALL ABORT
904   PRINT *,'UFBCPY - OUTPUT MESG IS NOT OPEN'
      CALL ABORT
905   PRINT *,'UFBCPY - IN/OUT I-NODE  MISMATCH'
      CALL ABORT
      END                                                               
