      subroutine writdx(lunit,lun,lundx)
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
      character*6   adn30
      character*1   moct(24000)
      dimension     mbay(5000)
      equivalence   (moct(1),mbay(1))

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check units, get a dx table, and start a dx message
c  ---------------------------------------------------

      if(lunit.eq.lundx) goto 900
      call readdx(lunit,lun,lundx)
      call dxmini(lun,mbay,mbyt,mby4,mbya,mbyb,mbyd)

      lda = ldxa(idxv+1)
      ldb = ldxb(idxv+1)
      ldd = ldxd(idxv+1)
      l30 = ld30(idxv+1)

c  copy table a contents to a bufr dx message
c  ------------------------------------------

      do i=1,ntba(lun)
      if(mbyt+lda.gt.maxdx) then
         call msgwrt(lunit,mbay,mbyt)
         call dxmini(lun,mbay,mbyt,mby4,mbya,mbyb,mbyd)
      endif
      call ipkm(moct(mby4),3,iupm(moct(mby4),24)+lda)
      call ipkm(moct(mbya),1,iupm(moct(mbya), 8)+  1)
      mbit = 8*(mbyb-1)
      call pkc(taba(i,lun),lda,mbay,mbit)
      call pkb(          0,  8,mbay,mbit)
      call pkb(          0,  8,mbay,mbit)
      mbyt = mbyt+lda
      mbyb = mbyb+lda
      mbyd = mbyd+lda
      enddo

c  copy table b contents to a bufr dx message
c  ------------------------------------------

      do i=1,ntbb(lun)
      if(mbyt+ldb.gt.maxdx) then
         call msgwrt(lunit,mbay,mbyt)
         call dxmini(lun,mbay,mbyt,mby4,mbya,mbyb,mbyd)
      endif
      call ipkm(moct(mby4),3,iupm(moct(mby4),24)+ldb)
      call ipkm(moct(mbyb),1,iupm(moct(mbyb), 8)+  1)
      mbit = 8*(mbyd-1)
      call pkc(tabb(i,lun),ldb,mbay,mbit)
      call pkb(          0,  8,mbay,mbit)
      mbyt = mbyt+ldb
      mbyd = mbyd+ldb
      enddo

c  copy table d contents to a bufr dx message
c  ------------------------------------------

      do i=1,ntbd(lun)
      nseq = iupm(tabd(i,lun)(ldd+1:ldd+1),8)
      lend = ldd+1 + l30*nseq
      if(mbyt+lend.gt.maxdx) then
         call msgwrt(lunit,mbay,mbyt)
         call dxmini(lun,mbay,mbyt,mby4,mbya,mbyb,mbyd)
      endif
      call ipkm(moct(mby4),3,iupm(moct(mby4),24)+lend)
      call ipkm(moct(mbyd),1,iupm(moct(mbyd), 8)+   1)
      mbit = 8*(mbyt-4)
      call pkc(tabd(i,lun),ldd,mbay,mbit)
      call pkb(       nseq,  8,mbay,mbit)
         do j=1,nseq
         jj  = ldd+2 + (j-1)*2
         idn = iupm(tabd(i,lun)(jj:jj),16)
         call pkc(adn30(idn,l30),l30,mbay,mbit)
         enddo
      mbyt = mbyt+lend
      enddo

c  write the unwritten message
c  ---------------------------

      call msgwrt(lunit,mbay,mbyt)

      return
900   call abt('writdx - inpur and output unit must not be the same')
      end
