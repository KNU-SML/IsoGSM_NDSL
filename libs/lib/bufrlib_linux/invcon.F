      function invcon(nc,lun,inv1,inv2)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      real*8 val
      real*4 vals

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  check the inventory interval
c  ----------------------------

      if(inv1.le.0 .or. inv1.gt.nval(lun)) goto 99
      if(inv2.le.0 .or. inv2.gt.nval(lun)) goto 99

c  find an occurance of node in the window meeting this condition
c  --------------------------------------------------------------

      do invcon=inv1,inv2
      if(inv(invcon,lun).eq.nodc(nc)) then
         if(kons(nc).eq.1 .and. val(invcon,lun).eq.vals(nc)) return
         if(kons(nc).eq.2 .and. val(invcon,lun).ne.vals(nc)) return
         if(kons(nc).eq.3 .and. val(invcon,lun).lt.vals(nc)) return
         if(kons(nc).eq.4 .and. val(invcon,lun).gt.vals(nc)) return
      endif
      enddo

99    invcon = 0
      return
      end
