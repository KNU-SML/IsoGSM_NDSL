      subroutine elemdx(card,lun)
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
      character*80  card
      character*24  unit
      character*11  refr
      character*8   nemo
      character*4   scal
      character*3   bitw
      character*1   sign,tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  capture the various elements characteristics
c  --------------------------------------------

      nemo = card( 3:10)
      scal = card(14:17)
      refr = card(21:31)
      bitw = card(35:37)
      unit = card(41:64)

c  find the element tag in table b
c  -------------------------------

      call nemtab(lun,nemo,idsn,tab,iele)
      if(tab.ne.'B') goto 900

c  left justify and store characteristics
c  --------------------------------------

      call jstchr(unit)
      tabb(iele,lun)(71:94) = unit

      call jstnum(scal,sign,iret)
      if(iret.ne.0) goto 901
      tabb(iele,lun)(95:95) = sign
      tabb(iele,lun)(96:98) = scal

      call jstnum(refr,sign,iret)
      if(iret.ne.0) goto 902
      tabb(iele,lun)( 99: 99) = sign
      tabb(iele,lun)(100:109) = refr

      call jstnum(bitw,sign,iret)
      if(iret.ne.0  ) goto 903
      if(sign.eq.'-') goto 903
      tabb(iele,lun)(110:112) = bitw

      return

c  error exit
c  -----------

900   call abt('elemdx - undefined element: '//card)
901   call abt('elemdx - bad scale value:   '//card)
902   call abt('elemdx - bad reference val: '//card)
903   call abt('elemdx - bad bit width:     '//card)
      end
