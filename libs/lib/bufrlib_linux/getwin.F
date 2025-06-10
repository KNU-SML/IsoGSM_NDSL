      subroutine getwin(node,lun,iwin,jwin)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      real*8 val

c----------------------------------------------------------------------
cfpp$ expand (lstrpc)
c----------------------------------------------------------------------

      irpc = lstrpc(node,lun)

      if(irpc.eq.0) then
         iwin = invwin(node,lun,jwin,nval(lun))
         if(iwin.eq.0 .and. jwin.gt.1) return
         iwin = 1
         jwin = nval(lun)
         return
      else
         iwin = invwin(irpc,lun,jwin,nval(lun))
         if(iwin.eq.0) then
            return
         elseif(val(iwin,lun).eq.0.) then
            iwin = 0
            return
         endif
      endif

      jwin = invwin(irpc,lun,iwin+1,nval(lun))
      if(jwin.eq.0) call abt('getwin - missing bracket')

      return
      end
