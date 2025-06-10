      subroutine rdtree(lun)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrbit/ nbit(inv_),mbit(inv_)
      common /scaleu/ sclu(-1000:1000)                                        

      character*10 tag
      character*8  cval
      character*3  typ
      dimension    ival(inv_)
      equivalence  (cval,rval)
c     real*8       val,rval
      real*8       val,rval,ups,ten
      data         ten/10./

c-----------------------------------------------------------------------
cfpp$ expand (upbb)
c-----------------------------------------------------------------------
      mps(node) = 2**(ibt(node))-1
c     ups(node) = (ival(n)+irf(node))*sclu(-isc(node))
      ups(node) = (ival(n)+irf(node))*ten**(-isc(node))
c-----------------------------------------------------------------------

c  cycle through a subset setting up the user array
c  ------------------------------------------------

      mbit(1) = ibit
      nbit(1) = 0
      call rcstpl(lun)

c  unpack a subset into the user array
c  -----------------------------------

      do n=1,nval(lun)
      call upbb(ival(n),nbit(n),mbit(n),lun)
      enddo

c  convert the unpacked integers to the proper types
c  -------------------------------------------------

      do n=1,nval(lun)
      node = inv(n,lun)
      if(itp(node).eq.1) then
         val(n,lun) = ival(n)
      elseif(itp(node).eq.2) then
         if(ival(n).lt.mps(node)) val(n,lun) = ups(node)
      endif
      enddo

c  special treatment for characters
c  --------------------------------

      do n=1,nval(lun)
      node = inv(n,lun)
      if(itp(node).eq.3) then
         cval = ' '
         call upc(cval,nbit(n)/8,mbay(1,lun),mbit(n))
         val(n,lun) = rval
      endif
      enddo

      return
      end
