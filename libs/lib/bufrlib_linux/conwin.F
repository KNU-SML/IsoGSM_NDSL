      subroutine conwin(lun,inc1,inc2,nbmp)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      real*8 val
      real*4 vals

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  special cases
c  -------------

      if(ncon.eq.0) then
         inc1 = 1
         inc2 = nval(lun)
         return
      endif

      if(inc1.gt.1 .and. kons(ncon).eq.5) then
         call nxtwin(lun,inc1,inc2)
         return
      endif

c  evaluate conditions to see if any more cases
c  --------------------------------------------

10    do nc=1,ncon
      if(kons(nc).eq.5) then
         inc1 = invwin(nodc(nc),lun,inc1,nval(lun))
         call usrtpl(lun,inc1-1,nbmp)
         call newwin(lun,inc1,inc2)
      else
15       call getwin(nodc(nc),lun,inc1,inc2)
         if(inc1.eq.0 .and. nc.eq.1) return
         if(inc1.eq.0              ) goto10
         icon = invcon(nc,lun,inc1,inc2)
         if(icon.eq.0) goto 15
      endif
      enddo

      return
      end
