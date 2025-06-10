      function lstrpc(node,lun)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*10 tag
      character*3  typ

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      if(node.lt.inode(lun) .or. node.gt.isc(inode(lun))) goto 900

      nod = node

c  find this or the previous rpc node
c  ----------------------------------

10    if(typ(nod).ne.'RPC') then
         nod = jmpb(nod)
         if(nod.ne.0) goto 10
      endif

      lstrpc = nod

      return
900   print*,inode(lun),':',node
      call abt('lstrpc - table node is out of bounds')
      end
