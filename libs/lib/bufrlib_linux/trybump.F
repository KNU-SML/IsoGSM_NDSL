      subroutine trybump(lunit,lun,usr,i1,i2,io,iret)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      real*8 usr(i1,i2),val
      real*4 vals

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  see if there is a drp group involved
c  ------------------------------------

      ndrp = lstjpb(nods(1),lun,'DRP')
      if(ndrp.le.0) return

c  if so, clean it out, bump it to i2, and try ufbrw again
c  -------------------------------------------------------

      invn = invwin(ndrp,lun,1,nval(lun))
      val(invn,lun) = 0
      jnvn = invn+1
      do while(nint(val(jnvn,lun)).gt.0)
         jnvn = jnvn+nint(val(jnvn,lun))
      enddo
      do knvn=1,nval(lun)-jnvn+1
      inv(invn+knvn,lun) = inv(jnvn+knvn-1,lun)
      val(invn+knvn,lun) = val(jnvn+knvn-1,lun)
      enddo
      nval(lun) = nval(lun)-(jnvn-invn-1)
      call usrtpl(lun,invn,i2)
      call ufbrw(lun,usr,i1,i2,io,iret)

      return
900   call abt('trybump - attempt to bump non-zero rep factor')
      end
