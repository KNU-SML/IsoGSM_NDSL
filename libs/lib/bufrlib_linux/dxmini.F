      subroutine dxmini(lun,mbay,mbyt,mb4,mba,mbb,mbd)
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
      dimension     mbay(5000)

      data nmcode/99/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  initialize the message
c  ----------------------

      mbit = 0
      do i=1,5000
      mbay(i) = 0
      enddo

      ih   = 0
      id   = 0
      im   = 0
      iy   = 0
      mtyp = 11
      msbt = idxv
      nsub = 1
      idxs = msbt+1
      ldxs = nxstr(idxs)

      nby0 = 8
      nby1 = 18
      nby2 = 0
      nby3 = 7 + nxstr(idxs) + 1
      nby4 = 7
      nby5 = 4
      mbyt = nby0+nby1+nby2+nby3+nby4+nby5

      if(mod(nby3,2).ne.0) goto 900

c  section 0
c  ---------

      call pkc('BUFR' ,  4 , mbay,mbit)
      call pkb(  mbyt , 24 , mbay,mbit)
      call pkb(     2 ,  8 , mbay,mbit)

c  section 1
c  ---------

      call pkb(  nby1 , 24 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)
      call pkb(nmcode , 16 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)
      call pkb(  mtyp ,  8 , mbay,mbit)
      call pkb(  msbt ,  8 , mbay,mbit)
      call pkb(     2 ,  8 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)
      call pkb(    iy ,  8 , mbay,mbit)
      call pkb(    im ,  8 , mbay,mbit)
      call pkb(    id ,  8 , mbay,mbit)
      call pkb(    ih ,  8 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)
      call pkb(     0 ,  8 , mbay,mbit)

c  section 3
c  ---------

      call pkb(       nby3 ,   24 , mbay,mbit)
      call pkb(          0 ,    8 , mbay,mbit)
      call pkb(          1 ,   16 , mbay,mbit)
      call pkb(       2**7 ,    8 , mbay,mbit)
c     call pkc(dxstr(idxs) , ldxs , mbay,mbit)
      do i=1,ldxs
      call pkb(iupm(dxstr(idxs)(i:i),8),8,mbay,mbit)
      enddo
      call pkb(          0 ,    8 , mbay,mbit)

c  section 4
c  ---------

      mb4 = mbit/8+1
      call pkb(nby4 , 24 , mbay,mbit)
      call pkb(   0 ,  8 , mbay,mbit)
      mba = mbit/8+1
      call pkb(   0 ,  8 , mbay,mbit)
      mbb = mbit/8+1
      call pkb(   0 ,  8 , mbay,mbit)
      mbd = mbit/8+1
      call pkb(   0 ,  8 , mbay,mbit)

      if(mbit/8+nby5.ne.mbyt) goto 901

      return
900   call abt('dxmini - uneven section 3')
901   call abt('dxmini - bytcnt is off   ')
      end
