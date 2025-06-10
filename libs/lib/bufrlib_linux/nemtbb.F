      subroutine nemtbb(lun,itab,unit,iscl,iref,ibit)
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
      character*8   nemo

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      if(itab.le.0 .or. itab.gt.ntbb(lun)) goto 900

c  pull out table b information
c  ----------------------------

      idn  = idnb(itab,lun)
      nemo = tabb(itab,lun)( 7:14)
      unit = tabb(itab,lun)(71:94)
      iscl = vals(tabb(itab,lun)( 95: 98))
      iref = vals(tabb(itab,lun)( 99:109))
      ibit = vals(tabb(itab,lun)(110:112))

c  check table b contents
c  ----------------------

      if(idn.lt.ifxy('000000')) goto 901
      if(idn.gt.ifxy('063255')) goto 901

      if(iscl.lt.-999 .or. iscl.gt.999) goto 902
      if(iref.le.-1e9 .or. iref.ge.1e9) goto 903
      if(ibit.le.   0 .or. ibit.gt. 64) goto 904

      if(unit.eq.'CCITT IA5' .and. mod(ibit,8).ne.0) goto 905

      return
900   call abt('nemtbb - itab not in table b'         )
901   call abt('nemtbb - bad descriptor value: '//nemo)
902   call abt('nemtbb - bad scale value     : '//nemo)
903   call abt('nemtbb - bad reference value : '//nemo)
904   call abt('nemtbb - bad bit width       : '//nemo)
905   call abt('nemtbb - bad char bit width  : '//nemo)
      end
