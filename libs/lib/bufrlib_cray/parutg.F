      SUBROUTINE PARUTG(LUN,IO,UTG,NOD,KON,VAL,*)                       
                                                                        
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
      COMMON /UTGPRM/ PICKEY                                            
                                                                        
      CHARACTER*20  UTG,ATAG                                            
      CHARACTER*10  TAG                                                 
      CHARACTER*3   TYP,ATYP,BTYP                                       
      CHARACTER*1   COND(5)                                             
      DIMENSION     BTYP(8),IOK(8)                                      
      LOGICAL       PICKEY                                              
                                                                        
      DATA PICKEY /.TRUE./                                              
      DATA NCHK   / 8/                                                  
      DATA BTYP   /'SUB','SEQ','REP','RPC','RPS','DRB','DRP','DRS'/     
      DATA IOK    /  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,   0 ,   0 /     
      DATA LTG    /20/                                                  
                                                                        
C---------------------------------------------------------------------- 
      COND(1) = '='                                                     
      COND(2) = '!'                                                     
      COND(3) = '<'                                                     
      COND(4) = '>'                                                     
      COND(5) = '^'                                                     
      NCOND   = 5                                                       
C---------------------------------------------------------------------- 
                                                                        
      ATAG  = ' '                                                       
      ATYP  = ' '                                                       
      KON   = 0                                                         
      NOD   = 0                                                         
      VAL   = 0                                                         
                                                                        
C  PARSE THE TAG                                                        
C  -------------                                                        
                                                                        
      DO I=1,LTG                                                        
      IF(UTG(I:I).EQ.' ') GOTO 1                                        
      DO J=1,NCOND                                                      
      IF(UTG(I:I).EQ.COND(J)) THEN                                      
         KON = J                                                        
         ICV = I+1                                                      
         GOTO 1                                                         
      ENDIF                                                             
      ENDDO                                                             
      ATAG(I:I) = UTG(I:I)                                              
      ENDDO                                                             
                                                                        
C  FIND THE TAG IN THE SUBSET TABLE                                     
C  --------------------------------                                     
                                                                        
1     INOD = INODE(LUN)                                                 
      DO NOD=INOD,ISC(INOD)                                             
      IF(ATAG.EQ.TAG(NOD)) GOTO 2                                       
      ENDDO                                                             
                                                                        
C     IF(KON.EQ.0 .AND. (IO.EQ.0.OR.ATAG.EQ.'NUL'.OR..NOT.PICKEY)) THEN 
      IF(KON.EQ.0) THEN                                                 
         NOD = 0                                                        
         RETURN                                                         
      ELSE                                                              
         RETURN 1                                                       
      ENDIF                                                             
                                                                        
C  CHECK FOR A VALID NODE TYP                                           
C  --------------------------                                           
                                                                        
2     IF(KON.EQ.5) THEN                                                 
         IF(TYP(NOD-1).NE.'DRP' .AND. TYP(NOD-1).NE.'DRS') GOTO 901     
      ELSE                                                              
         ATYP = TYP(NOD)                                                
         DO I=1,NCHK                                                    
         IF(ATYP.EQ.BTYP(I) .AND. IO.NE.IOK(I)) GOTO 902                
         ENDDO                                                          
      ENDIF                                                             
                                                                        
C  IF A COND NODE GET THE COND VALUE                                    
C  ---------------------------------                                    
                                                                        
      IF(KON.NE.0) THEN                                                 
         CALL STRNUM(UTG(ICV:LTG),NUM)                                  
         IF(NUM.LT.0) GOTO 903                                          
         VAL = NUM                                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
900   PRINT *,'PARUTG - NO VALID TAG FOUND IN              :',UTG  
      CALL ABORT
901   PRINT *,'PARUTG - BUMP NODE MUST BE TYPE RPC(DRP)    :',UTG  
      CALL ABORT
902   PRINT *,'PARUTG - ILLEGAL NODE TYPE:'//ATYP//       ':',UTG  
      CALL ABORT
903   PRINT *,'PARUTG - BAD OR MISSING COND VALUE IN       :',UTG  
      CALL ABORT
      END                                                               
