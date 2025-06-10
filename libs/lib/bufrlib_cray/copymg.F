      SUBROUTINE COPYMG(LUNIN,LUNOT)                                    
                                                                        
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      TASKCOMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      TASKCOMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10) 
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
                                                                        
      CHARACTER*10 TAG                                                  
      CHARACTER*8  SUBSET                                               
      CHARACTER*3  TYP                                                  
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  CHECK THE FILE STATUSES                                              
C  -----------------------                                              
                                                                        
      CALL STATUS(LUNIN,LIN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IL.GT.0) GOTO 901                                              
      IF(IM.LE.0) GOTO 902                                              
                                                                        
      CALL STATUS(LUNOT,LOT,IL,IM)                                      
      IF(IL.EQ.0) GOTO 903                                              
      IF(IL.LT.0) GOTO 904                                              
      IF(IM.NE.0) GOTO 905                                              
                                                                        
C  MAKE SURE BOTH FILES HAVE THE SAME TABLES                            
C  -----------------------------------------                            
                                                                        
      SUBSET = TAG(INODE(LIN))                                          
      CALL NEMTBA(LOT,SUBSET,MSGT,INOD)                                 
      IF(INODE(LIN).NE.INOD) GOTO 906                                   
                                                                        
C  EVERYTHING OKAY, COPY A MESSAGE                                      
C  -------------------------------                                      
                                                                        
      MBYM = IUPB(MBAY(1,LIN),5,24)                                     
      CALL MSGWRT(LUNOT,MBAY(1,LIN),MBYM)                               
                                                                        
C  SET THE MESSAGE CONTROL WORLDS FOR LUNOT                             
C  ----------------------------------------                             
                                                                        
      NMSG (LOT) = NMSG(LOT) + 1                                        
      NSUB (LOT) = MSUB(LIN)                                            
      IDATE(LOT) = IDATE(LIN)                                           
      INODE(LOT) = INODE(LIN)                                           
                                                                        
C  NORMAL EXIT                                                          
C  -----------                                                          
                                                                        
      RETURN                                                            
                                                                        
C  ERROR EXITS                                                          
C  -----------                                                          
                                                                        
900   PRINT *,'COPYMG - INPUT FILE IS CLOSED                    '
      CALL ABORT
901   PRINT *,'COPYMG - INPUT FILE IS OPEN FOR OUTPUT           '
      CALL ABORT
902   PRINT *,'COPYMG - NO INPUT FILE MESSAGE OPEN              '
      CALL ABORT
903   PRINT *,'COPYMG - OUTPUT FILE IS CLOSED                   '
      CALL ABORT
904   PRINT *,'COPYMG - OUTPUT FILE IS OPEN FOR OUTPUT          '
      CALL ABORT
905   PRINT *,'COPYMG - OUTPUT FILE MESSAGE OPEN                '
      CALL ABORT
906   PRINT *,'COPYMG - INPUT/OUTPUT FILES HAVE DIFFERENT TABLES'
      CALL ABORT
      END                                                               
