      subroutine chekstab(lun)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)
      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*24  unit
      character*8   nemo,nems(tbd_)
      character*1   tab
      dimension     irps(tbd_),knts(tbd_)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  there must be entries in tables a, b, and d
c  -------------------------------------------

      if(ntba(lun).eq.0) goto 900
      if(ntbb(lun).eq.0) goto 901
      if(ntbd(lun).eq.0) goto 902

c  make sure each table a entry defined as a sequence
c  --------------------------------------------------

      do i=1,ntba(lun)
      nemo = taba(i,lun)(4:11)
      call nemtab(lun,nemo,idn,tab,iret)
      if(tab.ne.'D') goto 903
      enddo

c  check table b contents
c  ----------------------

      do itab=1,ntbb(lun)
      call nemtbb(lun,itab,unit,iscl,iref,ibit)
      enddo

c  check table d contnets
c  ----------------------

      do itab=1,ntbd(lun)
      call nemtbd(lun,itab,nseq,nems,irps,knts)
      enddo

      return
900   call abt('chekstab - empty table a')
901   call abt('chekstab - empty table b')
902   call abt('chekstab - empty table d')
903   call abt('chekstab - no sequence defined for table a: '//nemo)
      end
