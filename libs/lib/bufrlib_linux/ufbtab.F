      subroutine ufbtab(lunit,tab,i1,i2,iret,str)
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
      common /scaleu/ sclu(-1000:1000)                                        

      character*(*) str
      character*10  tag,tgs(100)
      character*8   subset,cval
      character*3   typ
      dimension     tab(i1,i2),jbay(5000)
      equivalence   (cval,rval)
      logical       openit
      real*8        val,tab
      real*4 vals

      data maxtg /100/

c-----------------------------------------------------------------------
cfpp$ expand (upbb,usrtpl)
c-----------------------------------------------------------------------
      mps(node) = 2**(ibt(node))-1
      ups(node) = (ival+irf(node))*sclu(-isc(node))
c-----------------------------------------------------------------------

      iret = 0
      irec = 0
      isub = 0

      do j=1,i2
      do i=1,i1
      tab(i,j) = 10e10
      enddo
      enddo

c  see if we need to open a file
c  -----------------------------

      call status(lunit,lun,il,im)
      openit = il.eq.0

      if(openit) then
         call openbf(lunit,'IN',lunit)
      else
         call rewnbf(lunit,0)
      endif

c  check for special tags in string
c  --------------------------------

      call parseq(str,tgs,maxtg,ntg)
      do i=1,ntg
      if(tgs(i).eq.'IREC') irec = i
      if(tgs(i).eq.'ISUB') isub = i
      enddo

c  read a message and parse a string
c  ---------------------------------

10    call readmg(lunit,subset,idate,mret)
      if(mret.ne.0) goto 25
      call string(str,lun,i1,0)
      if(irec.gt.0) nods(irec) = 0
      if(isub.gt.0) nods(isub) = 0

c  see if there is another subset in the message
c  ---------------------------------------------

15    if(nsub(lun).eq.msub(lun)) goto 10
      if(iret+1.gt.i2) call abt('ufbtab - tab too small')
      iret = iret+1

      do i=1,nnod
      nods(i) = abs(nods(i))
      enddo

c  parse the string nodes from a subset
c  ------------------------------------

      mbit = mbyt(lun)*8 + 16
      nbit = 0
      n = 1
      call usrtpl(lun,n,n)
20    if(n+1.le.nval(lun)) then
         n = n+1
         node = inv(n,lun)
         mbit = mbit+nbit
         nbit = ibt(node)
         if(itp(node).eq.1) then
            call upbb(ival,nbit,mbit,lun)
            call usrtpl(lun,n,ival)
         endif
         do i=1,nnod
         if(nods(i).eq.node) then
            call upbb(ival,nbit,mbit,lun)
            if(itp(node).eq.1) then
               tab(i,iret) = ival
            elseif(itp(node).eq.2) then
               if(ival.lt.mps(node)) tab(i,iret) = ups(node)
            elseif(itp(node).eq.3) then
               cval = ' '
               call upc(cval,nbit/8,mbay(1,lun),mbit)
               tab(i,iret) = rval
            endif
            nods(i) = -nods(i)
            goto 20
         endif
         enddo
         do i=1,nnod
         if(nods(i).gt.0) goto 20
         enddo
      endif

c  update the subset pointers before next read
c  -------------------------------------------

      ibit = mbyt(lun)*8
      call upb(nbyt,16,mbay(1,lun),ibit)
      mbyt(lun) = mbyt(lun) + nbyt
      nsub(lun) = nsub(lun) + 1
      if(irec.gt.0) tab(irec,iret) = nmsg(lun)
      if(isub.gt.0) tab(isub,iret) = nsub(lun)
      goto 15

c  leave the file as it was before
c  -------------------------------

25    if(openit) then
         call closbf(lunit)
      else
         call rewnbf(lunit,1)
      endif

      return
      end
