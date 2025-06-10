      subroutine inctab(atag,atyp,node)
#include "ufbprm.h"
      integer tab_,seq_,tba_,tbb_,tbd_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*(*) atag,atyp
      character*10  tag
      character*3   typ

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      ntab = ntab+1
      if(ntab.gt.maxtab) call abt('inctab - too many entries')

      tag(ntab) = atag
      typ(ntab) = atyp
      node = ntab

      return
      end
