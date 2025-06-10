      subroutine usrtpl(lun,invn,nbmp)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      parameter (maxinv=inv_)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)

      character*10 tag
      character*3  typ
      dimension    itmp(100),vtmp(100)
      logical      drp,drs,drb,drx
      real*8       val,vtmp

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     print*,'usrtpl:',lun,':',invn,':',nbmp,':',tag(inode(lun))

      if(nbmp.le.0) return

      maxtmp = 100
      drp = .false.
      drs = .false.
      drx = .false.

c  set up a node expansion
c  -----------------------

      if(invn.eq.1) then
         nodi = inode(lun)
         inv(1,lun) = nodi
         nval(lun)  = 1
         if(nbmp.ne.1) goto 900
      elseif(invn.gt.0 .and. invn.le.nval(lun)) then
         nodi = inv(invn,lun)
         drp  = typ(nodi) .eq. 'DRP'
         drs  = typ(nodi) .eq. 'DRS'
         drb  = typ(nodi) .eq. 'DRB'
         drx  = drp .or. drs .or. drb
         ival = val(invn,lun)
         jval = 2**ibt(nodi)-1
         val(invn,lun) = ival+nbmp
         if(drb.and.nbmp.ne.1) goto 900
         if(.not.drx         ) goto 901
         if(ival.lt.0.       ) goto 902
         if(ival+nbmp.gt.jval) goto 903
      else
         goto 904
      endif

c  recall a pre-fab node expansion segment
c  ---------------------------------------

      newn = 0
      n1 = iseq(nodi,1)
      n2 = iseq(nodi,2)

      if(n1.eq.0          ) goto 905
      if(n2-n1+1.gt.maxtmp) goto 906

      do n=n1,n2
      newn = newn+1
      itmp(newn) = jseq(n)
      vtmp(newn) = vali(jseq(n))
      if(vtmp(newn).gt.10e9) vtmp(newn) = 10e10
      enddo

c  move old nodes - store new ones
c  -------------------------------

      if(nval(lun)+newn*nbmp.gt.maxinv) print*,'@:',nval(lun)+newn*nbmp
      if(nval(lun)+newn*nbmp.gt.maxinv) goto 907

cdir$ ivdep
      do j=nval(lun),invn+1,-1
      inv(j+newn*nbmp,lun) = inv(j,lun)
      val(j+newn*nbmp,lun) = val(j,lun)
      enddo

      if(drp.or.drs) vtmp(1) = newn
      knvn = invn

      do i=1,nbmp
      do j=1,newn
      knvn = knvn+1
      inv(knvn,lun) = itmp(j)
      val(knvn,lun) = vtmp(j)
      enddo
      enddo

c  reset pointers and counters
c  ---------------------------

      nval(lun) = nval(lun) + newn*nbmp

c     print*,tag(inv(invn,lun)),' ',newn,' ',nbmp,' ',nval(lun)
c     do i=1,newn
c     print*,tag(itmp(i))
c     enddo


      if(drx) then
         node = nodi
         invr = invn
4        node = jmpb(node)
         if(node.gt.0) then
            if(itp(node).eq.0) then
               do invr=invr-1,1,-1
               if(inv(invr,lun).eq.node) then
                  val(invr,lun) = val(invr,lun)+newn*nbmp
                  goto 4
               endif
               enddo
               goto 909
            else
               goto 4
            endif
         endif
      endif

c  normal exit
c  -----------

      return

c  error exits
c  -----------

900   call abt('usrtpl - nbmp <> 1 for        : '//tag(nodi))
901   call abt('usrtpl - node not sub,drp,drs : '//tag(nodi))
902   call abt('usrtpl - negative rep factor  : '//tag(nodi))
903   call abt('usrtpl - rep factor overflow  : '//tag(nodi))
904   call abt('usrtpl - inventory index out of bounds     ')
905   call abt('usrtpl - unset expansion seg  : '//tag(nodi))
906   call abt('usrtpl - temp array overflow  : '//tag(nodi))
907   call abt('usrtpl - inventory overflow   : '//tag(nodi))
908   call abt('usrtpl - tpl cache overflow   : '//tag(nodi))
909   call abt('usrtpl - bad backup strategy  : '//tag(nodi))
      end
