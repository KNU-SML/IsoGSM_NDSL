      subroutine nemtbd(lun,itab,nseq,nems,irps,knts)
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
      character*8   nemo,nems,nemt,nemf
      character*1   tab
      dimension     nems(tbd_),irps(tbd_),knts(tbd_)
      logical       rep

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      if(itab.le.0 .or. itab.gt.ntbd(lun)) goto 900

      rep  = .false.

c  clear the return values
c  -----------------------

      nseq = 0

      do i=1,tbd_
      nems(i) = ' '
      irps(i) = 0
      knts(i) = 0
      enddo

c  parse the table d entry
c  -----------------------

      nemo = tabd(itab,lun)(7:14)
      idsc = idnd(itab,lun)
      call uptdd(itab,lun,0,ndsc)

      if(idsc.lt.ifxy('300000')) goto 901
      if(idsc.gt.ifxy('363255')) goto 901
c     if(ndsc.le.0             ) goto 902

      do j=1,ndsc
      if(nseq+1.gt.tbd_) goto 903
      call uptdd(itab,lun,j,idsc)
      call numtab(lun,idsc,nemt,tab,iret)
      if(tab.eq.'R') then
         if(rep) goto 904
         rep = .true.
         if(iret.lt.0) then
            irps(nseq+1) = 1
            knts(nseq+1) = abs(iret)
         elseif(iret.gt.0) then
            irps(nseq+1) = iret
         endif
      elseif(tab.eq.'F') then
         if(.not.rep) goto 904
         irps(nseq+1) = iret
         rep = .false.
      elseif(tab.eq.'D') then
         rep = .false.
         nseq = nseq+1
         nems(nseq) = nemt
      elseif(tab.eq.'B') then
         rep = .false.
         nseq = nseq+1
         if(nemt(1:1).eq.'.') then
            call uptdd(itab,lun,j+1,idsc)
            call numtab(lun,idsc,nemf,tab,iret)
            call rsvfvm(nemt,nemf)
            if(tab.ne.'B') goto 906
         endif
         nems(nseq) = nemt
      else
         goto 905
      endif
      enddo

      return
900   call abt('nemtbd - itab not in table d   '                )
901   call abt('nemtbd - bad descriptor value: '          //nemo)
902   call abt('nemtbd - zero length sequence: '          //nemo)
903   call abt('nemtbd - too many descriptors in seq: '   //nemo)
904   call abt('nemtbd - replicator out of order in seq: '//nemo)
905   call abt('nemtbd - bad descriptor in sequence: '    //nemo)
906   call abt('nemtbd - following value not from tableb:'//nemf)
      end
