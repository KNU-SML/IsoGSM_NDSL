      subroutine ufbget(lunit,tab,i1,iret,str)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
 
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrbit/ nbit(inv_),mbit(inv_)
      common /scaleu/ sclu(-1000:1000)                                        
 
      character*(*) str
      character*10  tag,tgs(100)
      character*8   subset,bval
      character*3   typ
      dimension     tab(i1)
      equivalence   (ival,rval),(jval,bval)
      logical       openit
      real*8        val
      real*4 vals
 
      data maxtg /100/
 
c-----------------------------------------------------------------------
cfpp$ expand (upbb,usrtpl,invwin)
c-----------------------------------------------------------------------
      mps(node) = 2**(ibt(node))-1
      ups(node) = (ival+irf(node))*sclu(-isc(node))
c-----------------------------------------------------------------------
 
      iret = 0
      bval = ' '
 
c  make sure a file/message is open for input
c  ------------------------------------------
 
      call status(lunit,lun,il,im)
      if(il.ge.0) goto 900
      if(im.eq.0) goto 901
 
c  see if there is another subset in the message
c  ---------------------------------------------
 
      if(nsub(lun).eq.msub(lun)) then
         iret = -1
         return
      endif
 
c  parse the string
c  ----------------
 
      call string(str,lun,i1,0)
 
c  create the template for the zero expansion of this subset
c  ---------------------------------------------------------
 
      n = 1
      call usrtpl(lun,n,n)
 
      nbit(1) = 0
      mbit(1) = mbyt(lun)*8 + 16
 
      do n=2,nval(lun)
      node = inv(n,lun)
      nbit(n) = ibt(node)
      mbit(n) = mbit(n-1)+nbit(n-1)
      if(itp(node).eq.1 .or. node.eq.nods(nnod)) then
         nval(lun) = n
         goto 20
      endif
      enddo
20    continue
 
c  unpack only the nodes found in the string
c  -----------------------------------------
 
      do i=1,nnod
      node = nods(i)
      invn = invwin(node,lun,1,nval(lun))
      if(invn.gt.0) then
         call upbb(ival,nbit(invn),mbit(invn),lun)
         if(itp(node).eq.1) then
            tab(i) = ival
         elseif(itp(node).eq.2) then
            if(ival.lt.mps(node)) tab(i) = ups(node)
         elseif(itp(node).eq.3) then
            ival = ishft(ival,64-nbit(invn))
            ival = ior(ival,ishft(jval,-nbit(invn)))
            tab(i) = rval
         endif
      else
         tab(i) = 10e10
      endif
      enddo
 
      return
900   call abt('file not open for input')
901   call abt('no message open        ')
      end
