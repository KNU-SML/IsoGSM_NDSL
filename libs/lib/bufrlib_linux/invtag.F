      function invtag(node,lun,inv1,inv2)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(tab_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
 
      character*10 tag,tagn
      character*3  typ
      real*8       val
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      invtag = 0
      if(node.eq.0) return
      tagn = tag(node)
 
c  search between inv1 and inv2
c  ----------------------------
 
10    do invtag=inv1,inv2
      if(tag(inv(invtag,lun)).eq.tagn) return
      enddo
 
      invtag = 0
      return
      end
