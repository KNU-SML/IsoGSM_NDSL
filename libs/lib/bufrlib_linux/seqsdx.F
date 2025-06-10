      subroutine seqsdx(card,lun)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /reptab/ idnr(5,2),typs(5,2),reps(5,2),lens(5)

      character*80  card,seqs
      character*12  atag,tags(tbd_)
      character*8   nemo,nema,nemb
      character*3   typs
      character*1   reps,tab

      data maxtgs /tbd_/
      data maxtag /12/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  find the sequence tag in table d and parse the sequence string
c  --------------------------------------------------------------

      nemo = card( 3:10)
      seqs = card(14:78)

      call nemtab(lun,nemo,idn,tab,iseq)
      call parseq(seqs,tags,maxtgs,ntag)
      if(tab.ne.'D') goto 900
      if(ntag.eq.0 ) goto 900

      do n=1,ntag
      atag = tags(n)
      irep = 0

c  check for replicator
c  --------------------

      do i=1,5
      if(atag(1:1).eq.reps(i,1)) then
         do j=2,maxtag
         if(atag(j:j).eq.reps(i,2)) then
            if(j.eq.maxtag) goto 901
            call strnum(atag(j+1:maxtag),numr)
            if(i.eq.1 .and. numr.le.0  ) goto 901
            if(i.eq.1 .and. numr.gt.255) goto 901
            if(i.ne.1 .and. numr.ne.0  ) goto 901
            atag = atag(2:j-1)
            irep = i
            goto 1
         endif
         enddo
         goto 901
      endif
      enddo

c  check for valid tag
c  -------------------

1     if(nemock(atag).ne.0) goto 901
      call nemtab(lun,atag,idn,tab,iret)
      if(iret.gt.0) then
         if(tab.eq.'B' .and. irep.ne.0) goto 902
         if(atag(1:1).eq.'.') then
            nemb = tags(n+1)
            call numtab(lun,idn,nema,tab,itab)
            call nemtab(lun,nemb,jdn,tab,iret)
            call rsvfvm(nema,nemb)           
            if(nema.ne.atag) goto 903
            if(n.gt.ntag ) goto 905
            if(tab.ne.'B') goto 906
         endif
      else
         goto 903
      endif

c  write the descriptor string into tabd array
c  -------------------------------------------

10    if(irep.gt.0) call pktdd(iseq,lun,idnr(irep,1)+numr,iret)
      if(iret.lt.0) goto 904
      call pktdd(iseq,lun,idn,iret)
      if(iret.lt.0) goto 904

      enddo

      return

c  error exits
c  -----------

900   call abt('seqsdx - undefined sequence: '             //   nemo)
901   call abt('seqsdx - bad tag in sequence: '            //tags(n))
902   call abt('seqsdx - replicated elements not allowed:' //tags(n))
903   call abt('seqsdx - undefined tag: '                  //tags(n))
904   call abt('seqsdx - too many descriptors in string:'  //   nemo)
905   call abt('seqsdx - following-value last in string:'  //   nema)
906   call abt('seqsdx - following value not from tableb:' //   nemb)
      end
