      subroutine ufbcpy(lubin,lubot)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      real*8 val

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  check the file statuses and i-node
c  ----------------------------------

      call status(lubin,lui,il,im)
      if(il.eq.0) goto 900
      if(im.eq.0) goto 901
      if(inode(lui).ne.inv(1,lui)) goto 902

      call status(lubot,luo,il,im)
      if(il.eq.0) goto 903
      if(im.eq.0) goto 904
      if(inode(lui).ne.inode(luo)) goto 905

c  everything okay copy user array from lui to luo
c  -----------------------------------------------

      nval(luo) = nval(lui)

      do n=1,nval(lui)
      inv(n,luo) = inv(n,lui)
      val(n,luo) = val(n,lui)
      enddo

      return
900   call abt('ufbcpy - input  file is not open             ')
901   call abt('ufbcpy - input  mesg is not open             ')
902   call abt('ufbcpy - input  i-node  mismatch             ')
903   call abt('ufbcpy - output file is not open             ')
904   call abt('ufbcpy - output mesg is not open             ')
905   call abt('ufbcpy - in/out i-node  mismatch             ')
      end
