      SUBROUTINE INCTAB(ATAG,ATYP,NODE)
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)

      CHARACTER*(*) ATAG,ATYP
      CHARACTER*10  TAG
      CHARACTER*3   TYP

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      NTAB = NTAB+1
      IF(NTAB.GT.MAXTAB) CALL ABORT('INCTAB - TOO MANY ENTRIES')

      TAG(NTAB) = ATAG
      TYP(NTAB) = ATYP
      NODE = NTAB

      RETURN
      END
