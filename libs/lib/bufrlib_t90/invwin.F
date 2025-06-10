      FUNCTION INVWIN(NODE,LUN,INV1,INV2)                               
                                                                        
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
      REAL*8       VAL                                                  
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
      IF(NODE.EQ.0) RETURN                                              
                                                                        
C  SEARCH BETWEEN INV1 AND INV2                                         
C  ----------------------------                                         
                                                                        
10    DO INVWIN=INV1,INV2                                               
      IF(INV(INVWIN,LUN).EQ.NODE) RETURN                                
      ENDDO                                                             
                                                                        
      INVWIN = 0                                                        
      RETURN                                                            
      END                                                               
