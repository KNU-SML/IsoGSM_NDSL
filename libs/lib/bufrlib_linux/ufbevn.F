      subroutine ufbevn(lunit,usr,i1,i2,i3,iret,str)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      character*(*) str
      character*10  tag
      character*3   typ
      dimension     usr(i1,i2,i3),invn(255)
      real*8        val,usr
      real*4 vals

      data bmiss /10e10/

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  check the file status and i-node
c  --------------------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(im.eq.0) goto 901
      if(inode(lun).ne.inv(1,lun)) goto 902

c  parse the input string
c  ----------------------

      call string(str,lun,i1,0)

c  set initial values for returning arguments
c  ------------------------------------------

      do k=1,i3
      do j=1,i2
      do i=1,i1
      usr(i,j,k) = bmiss
      enddo
      enddo
      enddo

      iret = 0

c  loop over cond windows
c  ----------------------

      inc1 = 1
      inc2 = 1

1     call conwin(lun,inc1,inc2,i2)
      if(nnod.eq.0) then
         iret = i2
         return
      elseif(inc1.eq.0) then
         return
      else
         do i=1,nnod
         if(nods(i).gt.0) then
            ins2 = inc1
            call getwin(nods(i),lun,ins1,ins2)
            if(ins1.eq.0) return
            goto 2
         endif
         enddo
         ins1 = inc1
         ins2 = inc2
      endif

c  read push down stack data into 3d arrays
c  ----------------------------------------

2     iret = iret+1
      if(iret.le.i2) then
         do i=1,nnod
         if(nods(i).gt.0) then
            nnvn = nvnwin(nods(i),lun,ins1,ins2,invn,i3)
            do n=1,nnvn
            usr(i,iret,n) = val(invn(n),lun)
            enddo
         endif
         enddo
      endif

c  decide what to do next
c  ----------------------

      call nxtwin(lun,ins1,ins2)
      if(ins1.gt.0 .and. ins1.lt.inc2) goto 2
      if(ncon.gt.0) goto 1

      return
900   call abt('ufbevn - file is closed                     ')
901   call abt('ufbevn - no message open                    ')
902   call abt('ufbevn - i-node mismatch                    ')
      end
