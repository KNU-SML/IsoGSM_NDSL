      subroutine nemtab(lun,nemo,idn,tab,iret)
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
      character*8   nemt
      character*1   tab
      logical       folval
 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
 
      folval = nemo(1:1).eq.'.'
      iret = 0
      tab = ' '
 
c  look for nemo in table b
c  ------------------------
 
      do 1 i=1,ntbb(lun)
      nemt = tabb(i,lun)(7:14)
      if(nemt.eq.nemo) then
         idn  = idnb(i,lun)
         tab  = 'B'
         iret = i
         return
      elseif(folval.and.nemt(1:1).eq.'.') then
         do j=2,len(nemt)
         if(nemt(j:j).ne.'.' .and. nemt(j:j).ne.nemo(j:j)) goto 1
         enddo
         idn  = idnb(i,lun)
         tab  = 'B'
         iret = i
         return
      endif
1     enddo

c  dont look in table d for following value-mnemonics
c  ---------------------------------------------------

      if(folval) return
 
c  look in table d if we got this far
c  ----------------------------------
 
      do i=1,ntbd(lun)
      nemt = tabd(i,lun)(7:14)
      if(nemt.eq.nemo) then
         idn  = idnd(i,lun)
         tab  = 'D'
         iret = i
         return
      endif
      enddo
 
      return
      end
