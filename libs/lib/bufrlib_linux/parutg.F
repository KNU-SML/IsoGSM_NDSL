      subroutine parutg(lun,io,utg,nod,kon,val,*)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /utgprm/ pickey

      character*20  utg,atag
      character*10  tag
      character*3   typ,atyp,btyp
      character*1   cond(5)
      dimension     btyp(8),iok(8)
      logical       pickey

      data nchk   / 8/
      data btyp   /'SUB','SEQ','REP','RPC','RPS','DRB','DRP','DRS'/
      data iok    /  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,  -1 ,   0 ,   0 /
      data ltg    /20/

c----------------------------------------------------------------------
      pickey = .false.
      cond(1) = '='
      cond(2) = '!'
      cond(3) = '<'
      cond(4) = '>'
      cond(5) = '^'
      ncond   = 5
c----------------------------------------------------------------------

      atag  = ' '
      atyp  = ' '
      kon   = 0
      nod   = 0
      val   = 0

c  parse the tag
c  -------------

      do i=1,ltg
      if(utg(i:i).eq.' ') goto 1
      do j=1,ncond
      if(utg(i:i).eq.cond(j)) then
         kon = j
         icv = i+1
         goto 1
      endif
      enddo
      atag(i:i) = utg(i:i)
      enddo

c  find the tag in the subset table
c  --------------------------------

1     inod = inode(lun)
      do nod=inod,isc(inod)
      if(atag.eq.tag(nod)) goto 2
      enddo

      if(kon.eq.0 .and. (io.eq.0.or.atag.eq.'NUL'.or..not.pickey)) then
c     if(kon.eq.0) then
         nod = 0
         return
      else
         print*,'abort trying to write a non-existant mnemonic:'//atag
         return 1
      endif

c  check for a valid node typ
c  --------------------------

2     if(kon.eq.5) then
         if(typ(nod-1).ne.'DRP' .and. typ(nod-1).ne.'DRS') goto 901
      else
         atyp = typ(nod)
         do i=1,nchk
         if(atyp.eq.btyp(i) .and. io.ne.iok(i)) goto 902
         enddo
      endif

c  if a cond node get the cond value
c  ---------------------------------

      if(kon.ne.0) then
         call strnum(utg(icv:ltg),num)
         if(num.lt.0) goto 903
         val = num
      endif

      return
900   call abt('parutg - no valid tag found in              :'//utg)
901   call abt('parutg - bump node must be type rpc(DRP)    :'//utg)
902   call abt('parutg - illegal node type:'//atyp//       ':'//utg)
903   call abt('parutg - bad or missing cond value in       :'//utg)
      end
