      subroutine uptdd(id,lun,ient,iret)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)
      common /dxtab / maxdx,idxv,nxstr(10),ldxa(10),ldxb(10),ldxd(10),
     .                ld30(10),dxstr(10)

      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      ldd = ldxd(idxv+1)+1

c  check if ient is in bounds
c  --------------------------

      ndsc = iupm(tabd(id,lun)(ldd:ldd),8)

      if(ient.eq.0) then
         iret = ndsc
         return
      elseif(ient.lt.0 .or. ient.gt.ndsc) then
         call abt('uptdd - ient out of range')
      endif

c  return the descriptor indicated by ient
c  ---------------------------------------

      idsc = ldd+1 + (ient-1)*2
      iret = iupm(tabd(id,lun)(idsc:idsc),16)

      return
      end
