      SUBROUTINE RDTREE(LUN)                                            
                                                                        
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
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
      COMMON /SCALEU/ SCLU(-1000:1000)                                        
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  BVAL                                                 
      CHARACTER*3  TYP                                                  
      DIMENSION    IVAL(INV_),RVAL(INV_)                              
      EQUIVALENCE  (IVAL(1),RVAL(1)),(JVAL,BVAL)                        
      REAL*8       VAL,RVAL                                             
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND (UPBB)                                                     
C-----------------------------------------------------------------------
      MPS(NODE) = 2**(IBT(NODE))-1                                      
      UPS(NODE) = (IVAL(N)+IRF(NODE))*SCLU(-ISC(NODE))                 
C-----------------------------------------------------------------------
                                                                        
      BVAL = ' '                                                        
                                                                        
C  CYCLE THROUGH A SUBSET SETTING UP THE USER ARRAY                     
C  ------------------------------------------------                     
                                                                        
      MBIT(1) = IBIT                                                    
      NBIT(1) = 0                                                       
      CALL RCSTPL(LUN)                                                  
                                                                        
C  UNPACK A SUBSET INTO THE USER ARRAY                                  
C  -----------------------------------                                  
                                                                        
      DO N=1,NVAL(LUN)                                                  
      CALL UPBB(IVAL(N),NBIT(N),MBIT(N),LUN)                            
      ENDDO                                                             
                                                                        
C  CONVERT THE UNPACKED INTEGERS TO THE PROPER TYPES                    
C  -------------------------------------------------                    
                                                                        
      DO N=1,NVAL(LUN)                                                  
      NODE = INV(N,LUN)                                                 
      IF(ITP(NODE).EQ.1) THEN                                           
         VAL(N,LUN) = IVAL(N)                                           
      ELSEIF(ITP(NODE).EQ.2) THEN                                       
         IF(IVAL(N).LT.MPS(NODE)) VAL(N,LUN) = UPS(NODE)                
      ELSEIF(ITP(NODE).EQ.3) THEN                                       
         IVAL(N) = ISHFT(IVAL(N),64-NBIT(N))                            
         IVAL(N) = IOR(IVAL(N),ISHFT(JVAL,-NBIT(N)))                    
         VAL(N,LUN) = RVAL(N)                                           
      ENDIF                                                             
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
