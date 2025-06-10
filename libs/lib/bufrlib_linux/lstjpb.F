      function lstjpb(node,lun,jbtyp)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*(*) jbtyp
      character*10 tag
      character*3  typ

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  make sure we are all on the same page
c  -------------------------------------

      if(node.lt.inode(lun) .or. node.gt.isc(inode(lun))) then
         print*,inode(lun),':',node,':',tag(node)
         call abt('lstjpb - table node is out of bounds')
      endif

c  find this or the previous rpc node
c  ----------------------------------

      lstjpb = node

10    if(typ(lstjpb).ne.jbtyp) then
         lstjpb = jmpb(lstjpb)
         if(lstjpb.ne.0) goto 10
      endif

      return
      end
