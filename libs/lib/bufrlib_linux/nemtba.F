      subroutine nemtba(lun,nemo,mtyp,msbt,inod)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)

      character*(*) nemo
      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*20  nemt

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      nemt = nemo
      iret = 0

c  look for nemo in table a
c  ------------------------

      do i=1,ntba(lun)
      if(taba(i,lun)(4:11).eq.nemo) then
         mtyp = idna(i,lun,1)
         msbt = idna(i,lun,2)
         inod = mtab(i,lun)
         if(mtyp.lt.0 .or. mtyp.gt.255) goto 900
         if(msbt.lt.0 .or. msbt.gt.255) goto 901
         return
      endif
      enddo

      call abt('nemtba - cant find '//nemt)
900   call abt('nemtba - bad mtyp  '//nemt)
901   call abt('nemtba - bad msbt  '//nemt)
      end
