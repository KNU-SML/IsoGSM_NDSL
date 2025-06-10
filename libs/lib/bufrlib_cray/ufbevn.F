      SUBROUTINE UFBEVN(LUNIT,USR,I1,I2,I3,IRET,STR)                    
                                                                        
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
      TASKCOMMON /USRSTR/ NNOD,NCON,NODS(20),NODC(10),VALS(10),KONS(10) 
                                                                        
      CHARACTER*(*) STR                                                 
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP                                                 
      DIMENSION     USR(I1,I2,I3),INVN(255)                             
      REAL*8        VAL                                                 
                                                                        
      DATA BMISS /10E10/                                                
                                                                        
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
                                                                        
C  CHECK THE FILE STATUS AND I-NODE                                     
C  --------------------------------                                     
                                                                        
      CALL STATUS(LUNIT,LUN,IL,IM)                                      
      IF(IL.EQ.0) GOTO 900                                              
      IF(IM.EQ.0) GOTO 901                                              
      IF(INODE(LUN).NE.INV(1,LUN)) GOTO 902                             
                                                                        
C  PARSE THE INPUT STRING                                               
C  ----------------------                                               
                                                                        
      CALL STRING(STR,LUN,I1,0)                                         
                                                                        
C  SET INITIAL VALUES FOR RETURNING ARGUMENTS                           
C  ------------------------------------------                           
                                                                        
      USR = BMISS                                                       
      IRET = 0                                                          
                                                                        
C  LOOP OVER COND WINDOWS                                               
C  ----------------------                                               
                                                                        
      INC1 = 1                                                          
      INC2 = 1                                                          
                                                                        
1     CALL CONWIN(LUN,INC1,INC2,I2)                                     
      IF(NNOD.EQ.0) THEN                                                
         IRET = I2                                                      
         RETURN                                                         
      ELSEIF(INC1.EQ.0) THEN                                            
         RETURN                                                         
      ELSE                                                              
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            INS2 = INC1                                                 
            CALL GETWIN(NODS(I),LUN,INS1,INS2)                          
            GOTO 2                                                      
         ENDIF                                                          
         ENDDO                                                          
         INS1 = INC1                                                    
         INS2 = INC2                                                    
      ENDIF                                                             
                                                                        
C  READ PUSH DOWN STACK DATA INTO 3D ARRAYS                             
C  ----------------------------------------                             
                                                                        
2     IRET = IRET+1                                                     
      IF(IRET.LE.I2) THEN                                               
         DO I=1,NNOD                                                    
         IF(NODS(I).GT.0) THEN                                          
            NNVN = NVNWIN(NODS(I),LUN,INS1,INS2,INVN,I3)                
            DO N=1,NNVN                                                 
            USR(I,IRET,N) = VAL(INVN(N),LUN)                            
            ENDDO                                                       
         ENDIF                                                          
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  DECIDE WHAT TO DO NEXT                                               
C  ----------------------                                               
                                                                        
      CALL NXTWIN(LUN,INS1,INS2)                                        
      IF(INS1.GT.0 .AND. INS1.LT.INC2) GOTO 2                           
      IF(NCON.GT.0) GOTO 1                                              
                                                                        
      RETURN                                                            
900   PRINT *,'UFBINT - FILE IS CLOSED                     '
      CALL ABORT
901   PRINT *,'UFBINT - NO MESSAGE OPEN                    '
      CALL ABORT
902   PRINT *,'UFBINT - I-NODE MISMATCH                    '
      CALL ABORT
      END                                                               
