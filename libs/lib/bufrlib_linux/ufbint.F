      subroutine ufbint(lunin,usr,i1,i2,iret,str)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      character*(*) str
      dimension     usr(i1,i2)
      real*8        usr,val

c----------------------------------------------------------------------
cfpp$ expand (status,ufbrw)
c----------------------------------------------------------------------

c  check the file status and i-node
c  --------------------------------

      lunit = abs(lunin)
      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(im.eq.0) goto 901
      if(inode(lun).ne.inv(1,lun)) goto 902

      io = min(max(0,il),1)
      if(lunit.ne.lunin) io = 0

c  parse the input string - read/write values
c  ------------------------------------------

      call string(str,lun,i1,io)
      call ufbrw(lun,usr,i1,i2,io,iret)

      if(io.eq.1 .and. iret.ne.i2) then
         call trybump(lunit,lun,usr,i1,i2,io,iret)
         if(iret.ne.i2) print*,str
         if(iret.ne.i2) goto 903
      endif

      return
900   call abt('ufbint - file is closed                     ')
901   call abt('ufbint - no message open                    ')
902   call abt('ufbint - i-node mismatch                    ')
903   call abt('ufbint - incomplete write                   ')
      end
