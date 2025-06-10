      subroutine rdusdx(lundx,lun)
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
      character*8   nemo
      character*6   numb
      logical       nemtyp

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  initialize the dx-table partition and source file
c  -------------------------------------------------

      call dxinit(lun,1)
      rewind lundx

c  read user cards until there are no more
c  ---------------------------------------

1     read(lundx,'(a80)',end=100) card

c  reread if not a definition card
c  -------------------------------

      if(card(1: 1).eq.       '*') goto 1
      if(card(3:10).eq.'--------') goto 1
      if(card(3:10).eq.'        ') goto 1
      if(card(3:10).eq.'MNEMONIC') goto 1
      if(card(3:10).eq.'TABLE  D') goto 1
      if(card(3:10).eq.'TABLE  B') goto 1

c  parse a descriptor definition card
c  ----------------------------------

      if(card(12:12).eq.'|' .and. card(21:21).eq.'|') then

         nemo = card(3:10)
         numb = card(14:19)
         if(nemock(nemo).ne.0) goto 900
         if(numbck(numb).ne.0) goto 900

         if(numb(1:1).eq.'A') then
            n = ntba(lun)+1
            if(n.gt.ntba(0)) goto 901
            call nenuaa(nemo,numb(4:6),lun)
            taba(n,lun)( 1: 3) = numb(4:6)
            taba(n,lun)( 4:11) = nemo
            taba(n,lun)(13:67) = card(23:77)
            ntba(lun) = n

            read(nemo,'(2x,2i3)',err=10) mtyp,msbt
            idna(n,lun,1) = mtyp
            idna(n,lun,2) = msbt
            nemtyp = .true.
10          if(.not.nemtyp) then
               read(numb(4:6),'(i3)') idna(n,lun,1)
               idna(n,lun,2) = 0
            endif

            numb(1:1) = '3'
         endif

         if(numb(1:1).eq.'0') then
            n = ntbb(lun)+1
            if(n.gt.ntbb(0)) goto 902
            call nenubd(nemo,numb,lun)
            idnb(n,lun) = ifxy(numb)
            tabb(n,lun)( 1: 6) = numb
            tabb(n,lun)( 7:14) = nemo
            tabb(n,lun)(16:70) = card(23:77)
            ntbb(lun) = n
            goto 1
         endif

         if(numb(1:1).eq.'3') then
            n = ntbd(lun)+1
            if(n.gt.ntbd(0)) goto 903
            call nenubd(nemo,numb,lun)
            idnd(n,lun) = ifxy(numb)
            tabd(n,lun)( 1: 6) = numb
            tabd(n,lun)( 7:14) = nemo
            tabd(n,lun)(16:70) = card(23:77)
            ntbd(lun) = n
            goto 1
         endif

         goto 904
      endif

c  parse a sequence definition card
c  --------------------------------

      if(card(12:12).eq.'|' .and. card(19:19).ne.'|') then
         call seqsdx(card,lun)
         goto 1
      endif

c  parse an element definition card
c  --------------------------------

      if(card(12:12).eq.'|' .and. card(19:19).eq.'|') then
         call elemdx(card,lun)
         goto 1
      endif

c  cant figure out what kind of card it is
c  ----------------------------------------

      goto 905

c  normal exit
c  -----------

100   call makestab
      return

c  error exit
c  ----------

900   print*,card
      call abt('rdusdx - nemo or numb error             '//card)
901   call abt('rdusdx - too many table a entries       '//card)
902   call abt('rdusdx - too many table b entries       '//card)
903   call abt('rdusdx - too many table d entries       '//card)
904   call abt('rdusdx - bad descriptor number          '//card)
905   call abt('rdusdx - bad card format                '//card)
      end
