      SUBROUTINE WRTREE(LUN)                                            
                                                                        
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)                          
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
      DIMENSION    IVAL(INV_)                                          
      EQUIVALENCE  (RCHAR,ICHAR)                                        
      REAL*8       VAL                                                  
                                                                        
C-----------------------------------------------------------------------
      PKS(NODE) = VAL(N,LUN)*10.**ISC(NODE)-IRF(NODE)             
C-----------------------------------------------------------------------
                                                                        
C  CONVERT USER NUMBERS INTO SCALED INTEGERS                            
C  -----------------------------------------                            
                                                                        
      DO N=1,NVAL(LUN)                                                  
      NODE = INV(N,LUN)                                                 
      IF(ITP(NODE).EQ.1) THEN                                           
         IVAL(N) = VAL(N,LUN)                                           
      ELSEIF(ITP(NODE).EQ.2) THEN                                       
         IVAL(N) = NINT(PKS(NODE))
      ELSEIF(ITP(NODE).EQ.3) THEN                                       
         RCHAR   = VAL(N,LUN)                                           
         IVAL(N) = ISHFT(ICHAR,-(NBITW-IBT(NODE)))                      
      ENDIF                                                             
      ENDDO                                                             
                                                                        
C  PACK THE USER ARRAY INTO THE SUBSET BUFFER                           
C  ------------------------------------------                           
                                                                        
      IBIT = 16                                                         
                                                                        
      DO N=1,NVAL(LUN)                                                  
      NODE = INV(N,LUN)                                                 
      CALL PKB(IVAL(N),IBT(NODE),IBAY,IBIT)                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
