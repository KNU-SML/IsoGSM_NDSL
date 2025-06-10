      subroutine tabent(lun,nemo,tab,itab,irep,iknt,jum0)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)

      common /reptab/ idnr(5,2),typs(5,2),reps(5,2),lens(5)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*24 unit
      character*10 tag,rtag
      character*8  nemo
      character*3  typ,typs,typt
      character*1  reps,tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  make a jump/link table entry for a replicator
c  ---------------------------------------------

      if(irep.ne.0) then
         rtag = reps(irep,1)//nemo
         do i=1,10
         if(rtag(i:i).eq.' ') then
            rtag(i:i) = reps(irep,2)
            call inctab(rtag,typs(irep,1),node)
            jump(node) = node+1
            jmpb(node) = jum0
            link(node) = 0
            ibt (node) = lens(irep)
            irf (node) = 0
            isc (node) = 0
            if(irep.eq.1) irf(node) = iknt
            jum0 = node
            goto 1
         endif
         enddo
         goto 900
      endif

c  make an jump/link entry for an element or a sequence
c  ----------------------------------------------------

1     if(tab.eq.'B') then
         call nemtbb(lun,itab,unit,iscl,iref,ibit)
         if(unit.eq.'CCITT IA5') typt = 'CHR'
         if(unit.ne.'CCITT IA5') typt = 'NUM'
         call inctab(nemo,typt,node)
         jump(node) = 0
         jmpb(node) = jum0
         link(node) = 0
         ibt (node) = ibit
         irf (node) = iref
         isc (node) = iscl
      elseif(tab.eq.'D') then
         if(irep.eq.0) typt = 'SEQ'
         if(irep.ne.0) typt = typs(irep,2)
         call inctab(nemo,typt,node)
         jump(node) = node+1
         jmpb(node) = jum0
         link(node) = 0
         ibt (node) = 0
         irf (node) = 0
         isc (node) = 0
      else
         goto 901
      endif

      return
900   call abt('tabent - replicator error: '//rtag//':'//nemo)
901   call abt('tabent - undefined tag   : '           //nemo)
      end
