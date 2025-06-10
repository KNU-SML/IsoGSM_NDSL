      subroutine ufbrep(lunio,usr,i1,i2,iret,str)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
 
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
 
      character*(*) str
      dimension     usr(i1,i2)
      real*8        usr,val
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
c  check the file status and i-node
c  --------------------------------
 
      lunit = abs(lunio)
      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(im.eq.0) goto 901
      if(inode(lun).ne.inv(1,lun)) goto 902
      io = min(max(0,il),1)
      if(lunio.ne.lunit) io = 0
 
c  parse the input string - read/write values
c  ------------------------------------------
 
      call string(str,lun,i1,io)
      call ufbrp(lun,usr,i1,i2,io,iret)
 
      if(io.eq.1 .and. iret.ne.i2) then
         print*,str
         goto 903
      endif
 
      return
900   call abt('ufbrep - file is closed                     ')
901   call abt('ufbrep - no message open                    ')
902   call abt('ufbrep - i-node mismatch                    ')
903   call abt('ufbrep - incomplete write                   ')
      end
