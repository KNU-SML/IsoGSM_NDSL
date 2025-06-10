      subroutine dxinit(lun,ioi)
#include "ufbprm.h"
      integer tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /padesc/ ibct,ipd1,ipd2,ipd3,ipd4
      common /reptab/ idnr(5,2),typs(5,2),reps(5,2),lens(5)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)

      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*8   inib(6,5),inid(5)
      character*6   adn30
      character*3   typs
      character*1   reps

      data inib   /'------','BYTCNT  ','BYTES  ','+0','+0','16',
     .             '------','BITPAD  ','NONE   ','+0','+0','1 ',
     .             '031000','DRF1BIT ','NUMERIC','+0','+0','1 ',
     .             '031001','DRF8BIT ','NUMERIC','+0','+0','8 ',
     .             '031002','DRF16BIT','NUMERIC','+0','+0','16'/
      data ninib  /5/

      data inid   /'        ',
     .             'DRP16BIT',
     .             'DRP8BIT ',
     .             'DRPSTAK ',
     .             'DRP1BIT '/
      data ninid  /5/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  clear out a message control word partition
c  ------------------------------------------

      nmsg(lun)  = 0
      nsub(lun)  = 0
      msub(lun)  = 0
      inode(lun) = 0
      idate(lun) = 0

c  clear out a table partition
c  ---------------------------

      ntba(lun) = 0
      do i=1,ntba(0)
      taba(i,lun) = ' '
      mtab(i,lun) = 0
      enddo

      ntbb(lun) = 0
      do i=1,ntbb(0)
      tabb(i,lun) = ' '
      enddo

      ntbd(lun) = 0
      do i=1,ntbd(0)
      tabd(i,lun) = ' '
      call pktdd(i,lun,0,iret)
      enddo

      if(ioi.eq.0) return

c  initialize table with apriori table b and d entries
c  ---------------------------------------------------

      inib(1,1) = adn30(ibct,6)
      inib(1,2) = adn30(ipd4,6)

      do i=1,ninib
      ntbb(lun) = ntbb(lun)+1
      idnb(i,lun) = ifxy(inib(1,i))
      tabb(i,lun)(  1:  6) = inib(1,i)
      tabb(i,lun)(  7: 70) = inib(2,i)
      tabb(i,lun)( 71: 94) = inib(3,i)
      tabb(i,lun)( 95: 98) = inib(4,i)
      tabb(i,lun)( 99:109) = inib(5,i)
      tabb(i,lun)(110:112) = inib(6,i)
      enddo

      do i=2,ninid
      n = ntbd(lun)+1
      idnd(n,lun) = idnr(i,1)
      tabd(n,lun)(1: 6) = adn30(idnr(i,1),6)
      tabd(n,lun)(7:70) = inid(i)
      call pktdd(n,lun,idnr(1,1),iret)
      call pktdd(n,lun,idnr(i,2),iret)
      ntbd(lun) = n
      enddo

      return
      end
