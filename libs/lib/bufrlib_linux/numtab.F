      subroutine numtab(lun,idn,nemo,tab,iret)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /reptab/ idnr(5,2),typs(5,2),reps(5,2),lens(5)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)

      character*(*) nemo
      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*8   nemt
      character*3   typs
      character*1   reps,tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      nemo = ' '
      iret = 0
      tab = ' '

c  look for a replicator or a replicator factor
c  --------------------------------------------

      if(idn.ge.idnr(1,1) .and. idn.le.idnr(1,2)) then
         tab  = 'R'
         iret = -mod(idn,256)
         return
      endif

      do i=2,5
      if(idn.eq.idnr(i,1)) then
         tab  = 'R'
         iret = i
         return
      elseif(idn.eq.idnr(i,2)) then
         tab  = 'F'
         iret = i
         return
      endif
      enddo

c  look for idn in table d
c  -----------------------

      do i=1,ntbd(lun)
      if(idn.eq.idnd(i,lun)) then
         nemo = tabd(i,lun)(7:14)
         tab  = 'D'
         iret = i
         return
      endif
      enddo

c  look for idn in table b
c  -----------------------

      do i=1,ntbb(lun)
      if(idn.eq.idnb(i,lun)) then
         nemo = tabb(i,lun)(7:14)
         tab  = 'B'
         iret = i
         return
      endif
      enddo

      return
      end
