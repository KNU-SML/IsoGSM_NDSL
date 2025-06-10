      subroutine makestab
#include "ufbprm.h"
      integer tab_,seq_,tba_,tbb_,tbd_
      parameter(tba_=_tba_,tbb_=_tbb_,tbd_=_tbd_)
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)

      common /quiet/  iprt
      common /stbfr/  iolun(10),iomsg(10)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
      common /tababd/ ntba(0:10),ntbb(0:10),ntbd(0:10),mtab(tba_,10),
     .                idna(tba_,10,2),idnb(tbb_,10),idnd(tbd_,10),
     .                taba(tba_,10),tabb(tbb_,10),tabd(tbd_,10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)

      character*600 tabd
      character*128 tabb
      character*128 taba
      character*56  dxstr
      character*10  tag
      character*8   nemo
      character*3   typ
      dimension     lus(10)
      logical       expand,prttab
      real*8        val

      data prttab /.false./

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      prttab = iprt.ge.2

c  reset pointer table and string cache
c  ------------------------------------

      ntab = 0
      call strcln

c  figure out which units share tables
c  -----------------------------------

      do lun=1,10
      lus(lun) = 0
      if(iolun(lun).ne.0) then
         do lum=1,lun-1
         if(mtab(1,lun).eq.mtab(1,lum)) lus(lun) = lum
         enddo
      endif
      enddo

c  initialize jump-link tables with subsets/sequences/elements
c  ----------------------------------------------------------

      do lun=1,10

      if(iolun(lun).ne.0) then

c  reset any existing inventory pointers
c  -------------------------------------

         if(iomsg(lun).ne.0) then
            if(lus(lun).eq.0) inc = (ntab+1)-mtab(1,lun)
            if(lus(lun).ne.0) inc = mtab(1,lus(lun))-mtab(1,lun)
            do n=1,nval(lun)
            inv(n,lun) = inv(n,lun)+inc
            enddo
         endif

c  create new table entries if this unit doesnt share existing ones
c  -----------------------------------------------------------------

         if(lus(lun).eq.0) then
            call chekstab(lun)
            do itba=1,ntba(lun)
            inod = ntab+1
            nemo = taba(itba,lun)(4:11)
            call tabsub(lun,nemo)
            mtab(itba,lun) = inod
            isc(inod)      = ntab
c           do n1=inod,isc(inod)-1
c           do n2=n1+1,isc(inod)
c           if(tag(n1).eq.tag(n2)) goto 900
c           enddo
c           enddo
            enddo
         endif

      endif
      enddo

c  store types and initial values and counts
c  -----------------------------------------

      do node=1,ntab
      if(typ(node).eq.'SUB') then
         vali(node) = 0
         knti(node) = 1
         itp (node) = 0
      elseif(typ(node).eq.'SEQ') then
         vali(node) = 0
         knti(node) = 1
         itp (node) = 0
      elseif(typ(node).eq.'RPC') then
         vali(node) = 0
         knti(node) = 0
         itp (node) = 0
      elseif(typ(node).eq.'RPS') then
         vali(node) = 0
         knti(node) = 0
         itp (node) = 0
      elseif(typ(node).eq.'REP') then
         vali(node) = 10e10
         knti(node) = irf(node)
         itp (node) = 0
      elseif(typ(node).eq.'DRS') then
         vali(node) = 0
         knti(node) = 1
         itp (node) = 1
      elseif(typ(node).eq.'DRP') then
         vali(node) = 0
         knti(node) = 1
         itp (node) = 1
      elseif(typ(node).eq.'DRB') then
         vali(node) = 0
         knti(node) = 0
         itp (node) = 1
      elseif(typ(node).eq.'NUM') then
         vali(node) = 10e10
         knti(node) = 1
         itp (node) = 2
      elseif(typ(node).eq.'CHR') then
         vali(node) = 10e10
         knti(node) = 1
         itp (node) = 3
      else
         goto 901
      endif
      enddo

c  set up expansion segments for type 'SUB', 'DRP', and 'DRS' nodes
c  ----------------------------------------------------------------

      newn = 0

      do n=1,ntab
      iseq(n,1) = 0
      iseq(n,2) = 0
      expand = typ(n).eq.'SUB' .or. typ(n).eq.'DRP' .or. typ(n).eq.'DRS'
     .                         .or. typ(n).eq.'REP' .or. typ(n).eq.'DRB'
      if(expand) then
         iseq(n,1) = newn+1
         noda = n
         node = n+1
         do k=1,tab_
         knt(k) = 0
         enddo
         if(typ(noda).eq.'REP') knt(node) = knti(noda)
         if(typ(noda).ne.'REP') knt(node) = 1

1        newn = newn+1
         if(newn.gt.seq_) goto 902
         jseq(newn) = node
         knt(node) = max(knti(node),knt(node))
2        if(jump(node)*knt(node).gt.0) then
            node = jump(node)
            goto 1
         else if(link(node).gt.0) then
            node = link(node)
            goto 1
         else
            node = jmpb(node)
            if(node.eq.noda) goto 3
            if(node.eq.0   ) goto 903
            knt(node) = max(knt(node)-1,0)
            goto 2
         endif
3        iseq(n,2) = newn
      endif
      enddo

c  print the sequence tables
c  ------------------------

      if(prttab) then
         print*
         do i=1,ntab
         print99,i,
     .   tag(i),typ(i),jmpb(i),jump(i),link(i),ibt(i),irf(i),isc(i)
         enddo
         print*
99       format(i5,2x,a10,a5,6i8)
      endif

      return
900   call abt('makestab - dup in subset: '//tag(n1)//':'//nemo)
901   call abt('makestab - unknown type : '         //typ(node))
902   call abt('makestab - jseq overflow       : '  //tag(n   ))
903   call abt('makestab - failed to circulate : '  //tag(n   ))
      end
