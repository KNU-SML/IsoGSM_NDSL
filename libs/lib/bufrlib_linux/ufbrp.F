      subroutine ufbrp(lun,usr,i1,i2,io,iret)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
 
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(tab_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)
 
      character*10 tag
      character*3  typ
      real*8       usr(i1,i2),val
      real*4 vals
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      iret = 0
      ins1 = 0
      ins2 = 0
 
c  frame a section of the buffer - return when no frame
c  ----------------------------------------------------
 
1     if(ins1+1.gt.nval(lun)) return
      ins1 = invtag(nods(1),lun,ins1+1,nval(lun))
      if(ins1.eq.0) return
 
      ins2 = invtag(nods(1),lun,ins1+1,nval(lun))
      if(ins2.eq.0) ins2 = nval(lun)
      iret = iret+1
 
c  read user values
c  ----------------
 
      if(io.eq.0 .and. iret.le.i2) then
         do i=1,nnod
         usr(i,iret) = 10e10
         if(nods(i).gt.0) then
            invn = invtag(nods(i),lun,ins1,ins2)
            if(invn.gt.0) usr(i,iret) = val(invn,lun)
         endif
         enddo
      endif
 
c  write user values
c  -----------------
 
      if(io.eq.1 .and. iret.le.i2) then
         do i=1,nnod
         if(nods(i).gt.0) then
            invn = invtag(nods(i),lun,ins1,ins2)
            if(invn.gt.0) val(invn,lun) = usr(i,iret)
         endif
         enddo
      endif
 
c  go for next frame
c  -----------------
 
      goto 1
 
      end
