      subroutine pktdd(id,lun,idn,iret)
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

c  zero the counter if idn is zero
c  -------------------------------

      if(idn.eq.0) then
         call ipkm(tabd(id,lun)(ldd:ldd),1,0)
         iret = 0
         return
      endif

c  update the stored descriptor count for this table d entry
c  ---------------------------------------------------------

      nd = iupm(tabd(id,lun)(ldd:ldd),8)

      if(nd.lt.0 .or. nd.eq.tbd_) then
         iret = -1
         return
      else
         nd = nd+1
         call ipkm(tabd(id,lun)(ldd:ldd),1,nd)
         iret = nd
      endif

c  pack and store the descriptor
c  -----------------------------

      idm = ldd+1 + (nd-1)*2
      call ipkm(tabd(id,lun)(idm:idm),2,idn)

      return
      end
