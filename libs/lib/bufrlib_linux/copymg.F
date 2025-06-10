      subroutine copymg(lunin,lunot)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*10 tag
      character*8  subset
      character*3  typ

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check the file statuses
c  -----------------------

      call status(lunin,lin,il,im)
      if(il.eq.0) goto 900
      if(il.gt.0) goto 901
      if(im.le.0) goto 902

      call status(lunot,lot,il,im)
      if(il.eq.0) goto 903
      if(il.lt.0) goto 904
      if(im.ne.0) goto 905

c  make sure both files have the same tables
c  -----------------------------------------

      subset = tag(inode(lin))
      call nemtba(lot,subset,msgt,mstb,inod)
      if(inode(lin).ne.inod) goto 906

c  everything okay, copy a message
c  -------------------------------

      mbym = iupb(mbay(1,lin),5,24)
      call msgwrt(lunot,mbay(1,lin),mbym)

c  set the message control worlds for lunot
c  ----------------------------------------

      nmsg (lot) = nmsg(lot) + 1
      nsub (lot) = msub(lin)
      idate(lot) = idate(lin)
      inode(lot) = inode(lin)

c  normal exit
c  -----------

      return

c  error exits
c  -----------

900   call abt('copymg - input file is closed                    ')
901   call abt('copymg - input file is open for output           ')
902   call abt('copymg - no input file message open              ')
903   call abt('copymg - output file is closed                   ')
904   call abt('copymg - output file is open for output          ')
905   call abt('copymg - output file message open                ')
906   call abt('copymg - input/output files have different tables')
      end
