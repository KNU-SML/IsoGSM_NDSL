      subroutine nxtwin(lun,iwin,jwin)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      real*8 val

c----------------------------------------------------------------------
cfpp$ expand (lstrpc)
c----------------------------------------------------------------------

      if(jwin.eq.nval(lun)) then
         iwin = 0
         return
      endif

c  find the next sequential window
c  -------------------------------

      node = inv(iwin,lun)
      if(lstrpc(node,lun).ne.node) print*,'bad node=',node,iwin
      if(lstrpc(node,lun).ne.node) call abt('nxtwin - not rpc')
      if(val(jwin,lun).eq.0) then
         iwin = 0
      else
         iwin = jwin
         jwin = iwin+val(iwin,lun)
      endif

      return
      end
