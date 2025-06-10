      subroutine newwin(lun,iwin,jwin)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      real*8 val

c----------------------------------------------------------------------
cfpp$ expand (lstrpc)
c----------------------------------------------------------------------

      if(iwin.eq.1) then
         jwin = nval(lun)
         return
      endif

c  refind the jwin boundary from iwin
c  ----------------------------------

      node = inv(iwin,lun)
      if(lstrpc(node,lun).ne.node) call abt('newwin - not rpc')
      jwin = iwin+val(iwin,lun)

      return
      end
