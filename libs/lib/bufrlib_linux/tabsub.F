      subroutine tabsub(lun,nemo)
#include "ufbprm.h"
      integer tab_,seq_,tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)

      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*10 tag
      character*8  nemo,nems,nem
      character*3  typ
      character*1  tab
      dimension    nem(tbd_,10),irp(tbd_,10),krp(tbd_,10)
      dimension    drop(10),jmp0(10),nodl(10),ntag(10,2)
      logical      drop

      data maxlim /10/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check the mnemonic
c  ------------------

      call nemtab(lun,nemo,idn,tab,itab)
      if(tab.ne.'D') goto 900

c  store a subset node and jump/link the tree
c  ------------------------------------------

      call inctab(nemo,'SUB',node)
      jump(node) = node+1
      jmpb(node) = 0
      link(node) = 0
      ibt (node) = 0
      irf (node) = 0
      isc (node) = 0

      call nemtbd(lun,itab,nseq,nem(1,1),irp(1,1),krp(1,1))
      ntag(1,1) = 1
      ntag(1,2) = nseq
      jmp0(1)   = node
      limb      = 1

c  this loop resolves entities in a subset by emulating recursion
c  --------------------------------------------------------------

1     do n=ntag(limb,1),ntag(limb,2)

      ntag(limb,1) = n+1
      nodl(limb)   = ntab+1
      drop(limb)   = n.eq.ntag(limb,2)

      call nemtab(lun,nem(n,limb),idn,tab,itab)
      nems = nem(n,limb)
      irep = irp(n,limb)
      iknt = krp(n,limb)
      jum0 = jmp0(limb)
      call tabent(lun,nems,tab,itab,irep,iknt,jum0)

      if(tab.eq.'D') then
         limb = limb+1
         if(limb.gt.maxlim) goto 901
         call nemtbd(lun,itab,nseq,nem(1,limb),irp(1,limb),krp(1,limb))
         ntag(limb,1) = 1
         ntag(limb,2) = nseq
         jmp0(limb)   = ntab
         goto 1
      else if(drop(limb)) then
2        link(nodl(limb)) = 0
         limb = limb-1
         if(limb.eq.0 ) return
         if(drop(limb)) goto 2
         link(nodl(limb)) = ntab+1
         goto 1
      else
         link(nodl(limb)) = ntab+1
      endif

      enddo

      call abt('tabsub - should not get here               ')
900   call abt('tabsub - subset node not in table d: '//nemo)
901   call abt('tabsub - too many limbs                    ')
      end
