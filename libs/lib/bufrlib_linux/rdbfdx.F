      subroutine rdbfdx(lunit,lun)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)
      common /dxtab / maxdx,idxv,nxstr(10),ldxa(10),ldxb(10),ldxd(10),
     .                ld30(10),dxstr(10)

      character*600 tabd
      character*128 tabb,tabb1,tabb2
      character*128 taba
      character*56  dxstr
      character*50  dxcmp
      character*8   sec0,nemo,fmt30
      character*6   numb,cidn
      character*1   moct(24000)
      dimension     mbay(5000),ldxbd(10),ldxbe(10)
      equivalence   (mbay(1),moct(1))
      logical       nemtyp

      data ldxbd /38,70,8*0/
      data ldxbe /42,42,8*0/

c-----------------------------------------------------------------------
      ja(i) = ia+1+lda*(i-1)
      jb(i) = ib+1+ldb*(i-1)
c-----------------------------------------------------------------------

c  initialize the dx-table partition and source file
c  -------------------------------------------------

      call dxinit(lun,0)
      rewind lunit
      idx = 0

c  clear the buffer and read a message
c  -----------------------------------

1     do i=1,5000
      mbay(i) = 0
      enddo

      imsg = 8/nbytw+1
      read(lunit,err=900,end=2) sec0,(mbay(i),i=imsg,lmsg(sec0))
      idx = idx+1

c  get the section start octets and lengths
c  ----------------------------------------

2     i1 = 8
      l1 = iupm(moct(i1+1),24)
      i2 = i1+l1
      l2 = iupm(moct(i2+1),24)*iupm(moct(i1+8),1)
      i3 = i2+l2
      l3 = iupm(moct(i3+1),24)
      i4 = i3+l3
      l4 = iupm(moct(i4+1),24)

c  see if this is a bufr dx message - check for recognisable dx version
c  --------------------------------------------------------------------

      if(iupm(moct(i1+9),8).ne.11) then
c        backspace lunit

         rewind lunit
         do ndx=1,idx-1
         read(lunit,err=910,end=910) sec0,(mbay(i),i=imsg,lmsg(sec0))
         enddo

         call makestab
         return
      endif

      idxs = iupm(moct(i1+10),8)+1
      if(ldxa(idxs).eq.0) goto 902
      if(ldxb(idxs).eq.0) goto 902
      if(ldxd(idxs).eq.0) goto 902
      l30 = ld30(idxs)

      dxcmp = ' '
      call chrtrn(dxcmp,moct(i3+8),nxstr(idxs))
      if(dxcmp.ne.dxstr(idxs)) goto 902

c  section 4 - read definitions for tables a b and d
c  -------------------------------------------------

      lda  = ldxa (idxs)
      ldb  = ldxb (idxs)
      ldd  = ldxd (idxs)
      ldbd = ldxbd(idxs)
      ldbe = ldxbe(idxs)

      ia = i4+5
      la = iupm(moct(ia),8)
      ib = ja(la+1)
      lb = iupm(moct(ib),8)
      id = jb(lb+1)
      ld = iupm(moct(id),8)

c  table a - message type/subtype from the nemonic or the seq number
c  -----------------------------------------------------------------

      do i=1,la
      n = ntba(lun)+1
      if(n.gt.ntba(0)) goto 903
      call chrtrn(taba(n,lun),moct(ja(i)),lda)
      numb = taba(n,lun)(1:3)
      nemo = taba(n,lun)(4:11)
      call nenuaa(nemo,numb,lun)
      nemtyp = .false.
      ntba(lun) = n

      read(nemo,'(2x,2i3)',err=10) mtyp,msbt
      idna(n,lun,1) = mtyp
      idna(n,lun,2) = msbt
      nemtyp = .true.
10    if(.not.nemtyp) then
         read(numb(1:3),'(i3)') idna(n,lun,1)
         idna(n,lun,2) = 0
      endif
      enddo

c  table b
c  -------

      do i=1,lb
      n = ntbb(lun)+1
      if(n.gt.ntbb(0)) goto 904
      call chrtrna(tabb1,moct(jb(i)     ),ldbd)
      call chrtrna(tabb2,moct(jb(i)+ldbd),ldbe)
      tabb(n,lun) = tabb1(1:ldxbd(idxv+1))//tabb2(1:ldxbe(idxv+1))
      numb = tabb(n,lun)(1:6)
      nemo = tabb(n,lun)(7:14)
      call nenubd(nemo,numb,lun)
      idnb(n,lun) = ifxy(numb)
      ntbb(lun) = n
      enddo

c  table d
c  -------

      do i=1,ld
      n = ntbd(lun)+1
      if(n.gt.ntbd(0)) goto 905
      call chrtrna(tabd(n,lun),moct(id+1),ldd)
      numb = tabd(n,lun)(1:6)
      nemo = tabd(n,lun)(7:14)
      call nenubd(nemo,numb,lun)
      idnd(n,lun) = ifxy(numb)
      nd = iupm(moct(id+ldd+1),8)
      if(nd.gt.tbd_) goto 906
      do j=1,nd
      ndd = id+ldd+2 + (j-1)*l30
      call chrtrna(cidn,moct(ndd),l30)
      idn = idn30(cidn,l30)
      call pktdd(n,lun,idn,iret)
      if(iret.lt.0) goto 908
      enddo
      id = id+ldd+1 + nd*l30
      if(iupm(moct(id+1),8).eq.0) id = id+1
      ntbd(lun) = n
      enddo

c  goto read the next message
c  --------------------------

      goto 1

c  error exits
c  -----------

900   call abt('rdbfdx - i/o error reading dx message          ')
901   call abt('rdbfdx - eof >>>>> reading dx message          ')
902   call abt('rdbfdx - unexpected dx message type or contents')
903   call abt('rdbfdx - too many table a entries              ')
904   call abt('rdbfdx - too many table b entries              ')
905   call abt('rdbfdx - too many table d entries              ')
906   call abt('rdbfdx - too many descriptors in table d entry ')
907   call abt('rdbfdx - error reading idn seq from moct       ')
908   call abt('rdbfdx - bad return from pktdd                 ')
909   call abt('rdbfdx - desc count in tabd <> moct            ')
910   call abt('rdbfdx - err/eof positioning after dx messages ')
      end
