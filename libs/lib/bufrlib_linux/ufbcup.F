      subroutine ufbcup(lubin,lubot)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
 
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
 
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
 
      character*10 tag,tagi(inv_),tago
      character*3  typ
      dimension    nini(inv_)
      real*8       val
 
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
 
c  make a list of unique tags in input buffer
c  ------------------------------------------
 
      ntag = 0
 
      do 5 ni=1,nval(lui)
      nin = inv(ni,lui)
      if(itp(nin).ge.2) then
         do nv=1,ntag
         if(tagi(nv).eq.tag(nin)) goto 5
         enddo
         ntag = ntag+1
         nini(ntag) = ni
         tagi(ntag) = tag(nin)
      endif
5     enddo
 
      if(ntag.eq.0) goto 905
 
c  given a list make one copy of common elements to output buffer
c  --------------------------------------------------------------
 
      do 10 nv=1,ntag
      ni = nini(nv)
      do no=1,nval(luo)
      tago = tag(inv(no,luo))
      if(tagi(nv).eq.tago) then
         val(no,luo) = val(ni,lui)
         goto 10
      endif
      enddo
10    enddo
 
c  all exits here
c  --------------
 
      return
900   call abt('ufbcup - input  file is not open             ')
901   call abt('ufbcup - input  mesg is not open             ')
902   call abt('ufbcup - input  i-node  mismatch             ')
903   call abt('ufbcup - output file is not open             ')
904   call abt('ufbcup - output mesg is not open             ')
905   call abt('ufbcup - no tags in input buffer             ')
      end
