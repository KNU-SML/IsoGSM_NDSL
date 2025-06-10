      subroutine msgini(lun)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)

      common /padesc/ ibct,ipd1,ipd2,ipd3,ipd4
      common /msgptr/ nby0,nby1,nby2,nby3,nby4,nby5
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*10 tag
      character*8  subtag
      character*4  bufr,sevn
      character*3  typ
      character*1  tab

      data bufr/'BUFR'/
      data sevn/'7777'/
      data nmcd/   99 /

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  get the message tag and type, and break up the date
c  ---------------------------------------------------

      subtag = tag(inode(lun))
      call nemtba(lun,subtag,mtyp,msbt,inod)
      if(inode(lun).ne.inod) goto 900
      call nemtab(lun,subtag,isub,tab,iret)
      if(iret.eq.0) goto 901

      mear = mod(idate(lun)/1000000,100)
      mmon = mod(idate(lun)/10000  ,100)
      mday = mod(idate(lun)/100    ,100)
      mour = mod(idate(lun)/1      ,100)

c  initialize the message
c  ----------------------

      mbit = 0
      nby0 = 8
      nby1 = 18
      nby2 = 0
      nby3 = 20
      nby4 = 4
      nby5 = 4
      nbyt = nby0+nby1+nby2+nby3+nby4+nby5

c  section 0
c  ---------

      call pkc(bufr ,  4 , mbay(1,lun),mbit)
      call pkb(nbyt , 24 , mbay(1,lun),mbit)
      call pkb(   2 ,  8 , mbay(1,lun),mbit)

c  section 1
c  ---------

      call pkb(nby1 , 24 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(nmcd , 16 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(mtyp ,  8 , mbay(1,lun),mbit)
      call pkb(msbt ,  8 , mbay(1,lun),mbit)
      call pkb(   2 ,  8 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(mear ,  8 , mbay(1,lun),mbit)
      call pkb(mmon ,  8 , mbay(1,lun),mbit)
      call pkb(mday ,  8 , mbay(1,lun),mbit)
      call pkb(mour ,  8 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)

c  section 3
c  ---------

      call pkb(nby3 , 24 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)
      call pkb(   0 , 16 , mbay(1,lun),mbit)
      call pkb(2**7 ,  8 , mbay(1,lun),mbit)
      call pkb(ibct , 16 , mbay(1,lun),mbit)
      call pkb(isub , 16 , mbay(1,lun),mbit)
      call pkb(ipd1 , 16 , mbay(1,lun),mbit)
      call pkb(ipd2 , 16 , mbay(1,lun),mbit)
      call pkb(ipd3 , 16 , mbay(1,lun),mbit)
      call pkb(ipd4 , 16 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)

c  section 4
c  ---------

      call pkb(nby4 , 24 , mbay(1,lun),mbit)
      call pkb(   0 ,  8 , mbay(1,lun),mbit)

c  section 5
c  ---------

      call pkc(sevn ,  4 , mbay(1,lun),mbit)

c  double check initial message length
c  -----------------------------------

      if(mod(mbit,8).ne.0) goto 902
      if(mbit/8.ne.nbyt  ) goto 903

      nmsg(lun) = nmsg(lun)+1
      nsub(lun) = 0
      mbyt(lun) = nbyt

      return
900   call abt('msgini - inode mismatch                      ')
901   call abt('msgini - cant find subset idn: '//subtag      )
902   call abt('msgini - initial message out of byte boundary')
903   call abt('msgini - initial message count fails check   ')
      end
