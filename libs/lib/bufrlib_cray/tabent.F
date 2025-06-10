      SUBROUTINE TABENT(LUN,NEMO,TAB,ITAB,IREP,IKNT,JUM0)               
                                                                        
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)             
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
                                                                        
      CHARACTER*24 UNIT                                                 
      CHARACTER*10 TAG,RTAG                                             
      CHARACTER*8  NEMO                                                 
      CHARACTER*3  TYP,TYPS,TYPT                                        
      CHARACTER*1  REPS,TAB                                             
                                                                        
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                                                                        
C  MAKE A JUMP/LINK TABLE ENTRY FOR A REPLICATOR                        
C  ---------------------------------------------                        
                                                                        
      IF(IREP.NE.0) THEN                                                
         RTAG = REPS(IREP,1)//NEMO                                      
         DO I=1,10                                                      
         IF(RTAG(I:I).EQ.' ') THEN                                      
            RTAG(I:I) = REPS(IREP,2)                                    
            CALL INCTAB(RTAG,TYPS(IREP,1),NODE)                         
            JUMP(NODE) = NODE+1                                         
            JMPB(NODE) = JUM0                                           
            LINK(NODE) = 0                                              
            IBT (NODE) = LENS(IREP)                                     
            IRF (NODE) = 0                                              
            ISC (NODE) = 0                                              
            IF(IREP.EQ.1) IRF(NODE) = IKNT                              
            JUM0 = NODE                                                 
            GOTO 1                                                      
         ENDIF                                                          
         ENDDO                                                          
         GOTO 900                                                       
      ENDIF                                                             
                                                                        
C  MAKE AN JUMP/LINK ENTRY FOR AN ELEMENT OR A SEQUENCE                 
C  ----------------------------------------------------                 
                                                                        
1     IF(TAB.EQ.'B') THEN                                               
         CALL NEMTBB(LUN,ITAB,UNIT,ISCL,IREF,IBIT)                      
         IF(UNIT.EQ.'CCITT IA5') TYPT = 'CHR'                           
         IF(UNIT.NE.'CCITT IA5') TYPT = 'NUM'                           
         CALL INCTAB(NEMO,TYPT,NODE)                                    
         JUMP(NODE) = 0                                                 
         JMPB(NODE) = JUM0                                              
         LINK(NODE) = 0                                                 
         IBT (NODE) = IBIT                                              
         IRF (NODE) = IREF                                              
         ISC (NODE) = ISCL                                              
      ELSEIF(TAB.EQ.'D') THEN                                           
         IF(IREP.EQ.0) TYPT = 'SEQ'                                     
         IF(IREP.NE.0) TYPT = TYPS(IREP,2)                              
         CALL INCTAB(NEMO,TYPT,NODE)                                    
         JUMP(NODE) = NODE+1                                            
         JMPB(NODE) = JUM0                                              
         LINK(NODE) = 0                                                 
         IBT (NODE) = 0                                                 
         IRF (NODE) = 0                                                 
         ISC (NODE) = 0                                                 
      ELSE                                                              
         GOTO 901                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'TABENT - REPLICATOR ERROR: ',RTAG,':',NEMO
      CALL ABORT
901   PRINT *,'TABENT - UNDEFINED TAG   : '         ,NEMO
      CALL ABORT
      END                                                               
