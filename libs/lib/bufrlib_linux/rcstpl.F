      subroutine rcstpl(lun)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      parameter (maxtmp=1000)
      parameter (maxinv=inv_)
      parameter (maxrcr=100 )

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /usrbit/ nbit(inv_),mbit(inv_)

      character*10 tag
      character*3  typ
      dimension itmp(maxtmp,maxrcr),vtmp(maxtmp,maxrcr)
      dimension nbmp(2,maxrcr),newn(2,maxrcr)
      dimension knc(maxrcr),knx(maxrcr)
      real*8    val,vtmp

c-----------------------------------------------------------------------
cfpp$ expand (upbb)
c-----------------------------------------------------------------------

c  set the initial values for the template
c  ---------------------------------------

      inv(1,lun) = inode(lun)
      val(1,lun) = 0
      nbmp(1,1) = 1
      nbmp(2,1) = 1
      nodi = inode(lun)
      node = inode(lun)
      mbmp = 1
      knvn = 1
      nr   = 0

      do i=1,maxrcr
      knx(i) = 0
      enddo

c  set up the parametres for a level of recursion
c  ----------------------------------------------

10    continue

      nr = nr+1
      nbmp(1,nr) = 1
      nbmp(2,nr) = mbmp

      n1 = iseq(node,1)
      n2 = iseq(node,2)
      if(n1.eq.0          ) goto 905
      if(n2-n1+1.gt.maxtmp) goto 906
      newn(1,nr) = 1
      newn(2,nr) = n2-n1+1

      do n=1,newn(2,nr)
      nn = jseq(n+n1-1)
      itmp(n,nr) = nn
      vtmp(n,nr) = vali(nn)
      if(vtmp(n,nr).gt.10e9) vtmp(n,nr) = 10e10
      enddo

c  store nodes at some recursion level
c  -----------------------------------

20    do i=nbmp(1,nr),nbmp(2,nr)
      if(knx(nr).eq.0000) knx(nr) = knvn
      if(i.gt.nbmp(1,nr)) newn(1,nr) = 1
      do j=newn(1,nr),newn(2,nr)
      knvn = knvn+1
      node = itmp(j,nr)
      inv(knvn,lun) = node
      val(knvn,lun) = vtmp(j,nr)
      mbit(knvn) = mbit(knvn-1)+nbit(knvn-1)
      nbit(knvn) = ibt(node)
      if(itp(node).eq.1) then
         call upbb(mbmp,nbit(knvn),mbit(knvn),lun)
         newn(1,nr) = j+1
         nbmp(1,nr) = i
         goto 10
      endif
      enddo
      new = knvn-knx(nr)
      val(knx(nr)+1,lun) = val(knx(nr)+1,lun) + new
      knx(nr) = 0
      enddo

c  continue at one recusion level back
c  -----------------------------------

      if(nr-1.ne.0) then
         nr = nr-1
         goto 20
      endif

c  finally store the length of the subset template
c  -----------------------------------------------

      nval(lun) = knvn

c  normal exit
c  -----------

      return

c  error exits
c  -----------

900   call abt('rcstpl - nbmp <> 1 for        : '//tag(nodi))
901   call abt('rcstpl - node not sub,drp,drs : '//tag(nodi))
902   call abt('rcstpl - negative rep factor  : '//tag(nodi))
903   call abt('rcstpl - rep factor overflow  : '//tag(nodi))
904   call abt('rcstpl - inventory index out of bounds     ')
905   call abt('rcstpl - unset expansion seg  : '//tag(nodi))
906   call abt('rcstpl - temp array overflow  : '//tag(nodi))
907   call abt('rcstpl - inventory overflow   : '//tag(nodi))
908   call abt('rcstpl - tpl cache overflow   : '//tag(nodi))
909   call abt('rcstpl - bad backup strategy  : '//tag(nodi))
      end
