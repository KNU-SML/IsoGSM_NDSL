      subroutine ufbrw(lun,usr,i1,i2,io,iret)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      character*10 tag
      character*3  typ
      real*8       usr(i1,i2),val
      real*4 vals

c----------------------------------------------------------------------
cfpp$ expand (conwin,drstpl,getwin,invwin,lstrps,newin,nxtwin)
c----------------------------------------------------------------------

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

c  loop over store nodes
c  ---------------------

2     iret = iret+1

c     print*,'ufbrw:',iret,':',ins1,':',ins2,':',inc1,':',inc2
c     print'(5a10)',(tag(inv(i,lun)),i=ins1,ins2)

c  write user values
c  -----------------

      if(io.eq.1 .and. iret.le.i2) then
         do i=1,nnod
         if(nods(i).gt.0) then
            if(usr(i,iret).ne.10e10) then
               invn = invwin(nods(i),lun,ins1,ins2)
               if(invn.eq.0) then
                  call drstpl(nods(i),lun,ins1,ins2,invn)
                  if(invn.eq.0) then
                     iret = 0
                     return
                  endif
                  call newwin(lun,inc1,inc2)
                  val(invn,lun) = usr(i,iret)
               elseif(lstrps(nods(i),lun).eq.0) then
                  val(invn,lun) = usr(i,iret)
               elseif(val(invn,lun).eq.10e10) then
                  val(invn,lun) = usr(i,iret)
               else
                  call drstpl(nods(i),lun,ins1,ins2,invn)
                  if(invn.eq.0) then
                     iret = 0
                     return
                  endif
                  call newwin(lun,inc1,inc2)
                  val(invn,lun) = usr(i,iret)
               endif
            endif
         endif
         enddo
      endif

c  read user values
c  ----------------

      if(io.eq.0 .and. iret.le.i2) then
         do i=1,nnod
         usr(i,iret) = 10e10
         if(nods(i).gt.0) then
            invn = invwin(nods(i),lun,ins1,ins2)
            if(invn.gt.0) usr(i,iret) = val(invn,lun)
         endif
         enddo
      endif

c  decide what to do next
c  ----------------------

      if(io.eq.1.and.iret.eq.i2) return
      call nxtwin(lun,ins1,ins2)
      if(ins1.gt.0 .and. ins1.lt.inc2) goto 2
      if(ncon.gt.0) goto 1

      return
      end
