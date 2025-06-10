      subroutine nenuck(nemo,numb,lun)
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
      character*8   nemo
      character*6   numb

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check table a
c  -------------

      entry nenuaa(nemo,numb,lun)

      do n=1,ntba(lun)
      if(numb(2:4).eq.taba(n,lun)(1: 3)) goto 900
      if(nemo     .eq.taba(n,lun)(4:11)) goto 900
      enddo

      return

c  check table b and d
c  -------------------

      entry nenubd(nemo,numb,lun)

      do n=1,ntbb(lun)
      if(numb.eq.tabb(n,lun)(1: 6)) goto 900
      if(nemo.eq.tabb(n,lun)(7:14)) goto 900
      enddo

      do n=1,ntbd(lun)
      if(numb.eq.tabd(n,lun)(1: 6)) goto 900
      if(nemo.eq.tabd(n,lun)(7:14)) goto 900
      enddo

      return

c  error exit
c  ----------

900   call abt('nenuck - duplicate nem/num '//nemo//' '//numb)
      end
