      subroutine wrtree(lun)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /scaleu/ sclu(-1000:1000)                                        

      character*10 tag
      character*8  cval
      character*3  typ
      dimension    ival(inv_)
      equivalence  (cval,rval)
      real*8       val,rval,pks,ten
      data ten /10./

c-----------------------------------------------------------------------
c     pks(node) = val(n,lun)*sclu(isc(node))-irf(node)
      pks(node) = val(n,lun)*ten**(isc(node))-irf(node)
c-----------------------------------------------------------------------

c  convert user numbers into scaled integers
c  -----------------------------------------

      do n=1,nval(lun)
      node = inv(n,lun)
      if(itp(node).eq.1) then
         ival(n) = val(n,lun)
      elseif(itp(node).eq.2) then
         ival(n) = nint(pks(node))
      endif
      enddo

c  pack the user array into the subset buffer
c  ------------------------------------------

      ibit = 16

      do n=1,nval(lun)
      node = inv(n,lun)
      if(itp(node).lt.3) then
         call pkb(ival(n),ibt(node),ibay,ibit)
      else
         rval = val(n,lun)
         call pkc(cval,ibt(node)/8,ibay,ibit)
      endif
      enddo

      return
      end
