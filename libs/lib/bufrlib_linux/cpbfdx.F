      subroutine cpbfdx(lud,lun)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)

      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  initialize the dx-table partition
c  ---------------------------------

      call dxinit(lun,0)

c  copy one table partition to another
c  -----------------------------------

      inode(lun) = inode(lud)

      ntba(lun) = ntba(lud)
      ntbb(lun) = ntbb(lud)
      ntbd(lun) = ntbd(lud)

      do i=1,ntba(lud)
      idna(i,lun,1) = idna(i,lud,1)
      idna(i,lun,2) = idna(i,lud,2)
      taba(i,lun) = taba(i,lud)
      mtab(i,lun) = mtab(i,lud)
      enddo

      do i=1,ntbb(lud)
      idnb(i,lun) = idnb(i,lud)
      tabb(i,lun) = tabb(i,lud)
      enddo

      do i=1,ntbd(lud)
      idnd(i,lun) = idnd(i,lud)
      tabd(i,lun) = tabd(i,lud)
      enddo

      return
      end
