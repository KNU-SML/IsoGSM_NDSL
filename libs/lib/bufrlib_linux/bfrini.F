      subroutine bfrini
#include "ufbprm.h"
      integer tab_,seq_,tba_,tbb_,tbd_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)

      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /padesc/ ibct,ipd1,ipd2,ipd3,ipd4
      common /reptab/ idnr(5,2),typs(5,2),reps(5,2),lens(5)
      common /scaleu/ sclu(-1000:1000)                                        
      common /stbfr / iolun(10),iomsg(10)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)
      common /dxtab / maxdx,idxv,nxstr(10),ldxa(10),ldxb(10),ldxd(10),
     .                ld30(10),dxstr(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /bufrmg/ msglen,msgtxt(5000)


      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*10  tag
      character*6   adsn(5,2),dndx(25,10)
      character*3   typx(5,2),typs,typ
      character*1   repx(5,2),reps
      dimension     ndndx(10),nldxa(10),nldxb(10),nldxd(10),nld30(10)
      dimension     lenx(5)

      data adsn   / '101000','360001','360002','360003','360004' ,
     .              '101255','031002','031001','031001','031000' /
      data typx   /    'REP',   'DRP',   'DRP',   'DRS' ,  'DRB' ,
     .                 'SEQ',   'RPC',   'RPC',   'RPS' ,  'SEQ' /
      data repx   /      '"',     '(',     '{',     '[' ,    '<' ,
     .                   '"',     ')',     '}',     ']' ,    '>' /
      data lenx   /       0 ,     16 ,      8 ,      8  ,     1  /

      data (dndx(i,1),i=1,25)/
     .'102000','031001','000001','000002',
     .'110000','031001','000010','000011','000012','000013','000015',
     .                  '000016','000017','000018','000019','000020',
     .'107000','031001','000010','000011','000012','000013','101000',
     .                  '031001','000030'/

      data (dndx(i,2),i=1,15)/
     .'103000','031001','000001','000002','000003',
     .'101000','031001','300004',
     .'105000','031001','300003','205064','101000','031001','000030'/

      data ndndx /  25 ,  15 , 8*0 /
      data nldxa /  35 ,  67 , 8*0 /
      data nldxb /  80 , 112 , 8*0 /
      data nldxd /  38 ,  70 , 8*0 /
      data nld30 /   5 ,   6 , 8*0 /

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  initialize /bitbuf/
c  -------------------

c     maxbyt = 9970
      maxbyt=19950

c  initialize /padesc/
c  -------------------

      ibct = ifxy('063000')
      ipd1 = ifxy('102000')
      ipd2 = ifxy('031001')
      ipd3 = ifxy('206001')
      ipd4 = ifxy('063255')

c  initialize /stbfr/
c  ------------------

      do i=1,10
      iolun(i) = 0
      iomsg(i) = 0
      enddo

c  initialize /reptab/
c  -------------------

      do i=1,5
      lens(i) = lenx(i)
      do j=1,2
      idnr(i,j) = ifxy(adsn(i,j))
      typs(i,j) = typx(i,j)
      reps(i,j) = repx(i,j)
      enddo
      enddo

c  initialize /tababd/
c  -------------------

      ntba(0) = tba_
      ntbb(0) = tbb_
      ntbd(0) = tbd_

c  initialize /dxtab/
c  ------------------

      maxdx = maxbyt
      idxv  = 1

      do j=1,10
      ldxa(j)  = nldxa(j)
      ldxb(j)  = nldxb(j)
      ldxd(j)  = nldxd(j)
      ld30(j)  = nld30(j)
      dxstr(j) = '      '
      nxstr(j) = ndndx(j)*2
      do i=1,ndndx(j)
      i1 = i*2-1
      call ipkm(dxstr(j)(i1:i1),2,ifxy(dndx(i,j)))
      enddo
      enddo

c  initialize /tables/
c  -------------------

      maxtab = tab_

c  initialize /bufrmg/
c  -------------------

      msglen = 0
                                                                        
c  initialize /scaleu/                                                  
c  -------------------                                                  
                                                                        
c     sclu(-1) = .1
c     sclu( 0) =  1
c     sclu( 1) = 10
c     do i=2,1000
c     sclu(-i) = sclu(1-i)*sclu(-1)
c     sclu( i) = sclu(i-1)*sclu( 1)
c     enddo                            

      return
      end
