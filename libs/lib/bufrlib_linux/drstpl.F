      subroutine drstpl(inod,lun,inv1,inv2,invn)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      character*10 tag
      character*3  typ
      real*8       val

c-----------------------------------------------------------------------
cfpp$ expand (invwin,usrtpl,newwin)
c-----------------------------------------------------------------------

1     node = inod
2     node = jmpb(node)
      if(node.eq.0) return
      if(typ(node).eq.'DRS' .or. typ(node).eq.'DRB') then
         invn = invwin(node,lun,inv1,inv2)
         if(invn.gt.0) then
            call usrtpl(lun,invn,1)
            call newwin(lun,inv1,inv2)
            invn = invwin(inod,lun,invn,inv2)
            if(invn.gt.0) return
            goto 1
         endif
      endif
      goto 2
900   call abt('drstpl - cant find node:'//tag(inod))
c     print*,'drstpl:',tag(inod),':',tag(node),inv1,inv2
c     print'(5a10)',(tag(inv(i,lun)),i=inv1,inv2)
      end
