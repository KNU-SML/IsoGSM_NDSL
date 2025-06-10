      subroutine ufbdmp(lunit,luprt)
#include "ufbprm.h"
      integer tab_,seq_
      parameter(inv_=_inv_,tab_=_tab_,seq_=_seq_)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /tables/ maxtab,ntab,tag(tab_),typ(tab_),knt(tab_),
     .                jump(tab_),link(tab_),jmpb(tab_),
     .                ibt(tab_),irf(tab_),isc(tab_),
     .                itp(tab_),vali(tab_),knti(tab_),
     .                iseq(tab_,2),jseq(seq_)
      common /usrint/ nval(10),inv(inv_,10),val(inv_,10)
 
      character*10 tag,tg
      character*8  vc
      character*3  typ,tp
      character*1  you
      equivalence  (vl,vc)
      real*8       val,vl
 
      data bmiss /10e10/
      data xmsg  /'MISSING'/
 
c----------------------------------------------------------------------
cffp$ expand (status)
c----------------------------------------------------------------------
 
      if(luprt.eq.0) luout = 6
      if(luprt.ne.0) luout = luprt
 
c  check the file status and i-node
c  --------------------------------
 
      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(im.eq.0) goto 901
      if(inode(lun).ne.inv(1,lun)) goto 902
 
c  dump the contents of common /usrint/ for unit lunit
c  ---------------------------------------------------
 
      do nv=1,nval(lun)
      if(luprt.eq.0 .and. mod(nv,20).eq.0) then
         print*,'(MORE)'
         read(5,'(a1)') you
         if(you.eq.'q') return
      endif
      nd = inv (nv,lun)
      vl = val (nv,lun)
      tg = tag (nd)
      tp = typ (nd)
      it = itp (nd)
      ib = ibt (nd)
      is = isc (nd)
      ir = irf (nd)
      jp = jump(nd)
      lk = link(nd)
      jb = jmpb(nd)
      rj = rjust(tg)
      if(tp.ne.'CHR' .and. nint(vl).lt.bmiss) then
         write(luout,1) nv,tp,it,tg,vl,ib,is,ir,nd,jp,lk,jb
      else
         if(tp.ne.'CHR') vl = xmsg
         if(tp.eq.'CHR') vc = vc(1:ib/8)
         rj = rjust(vc)
         write(luout,2) nv,tp,it,tg,vc,ib,is,ir,nd,jp,lk,jb
      endif
      enddo
 
1     format(i5,1x,a3,'-',i1,1x,a10,1x,f10.1,7(1x,i5))
2     format(i5,1x,a3,'-',i1,1x,a10,1x,a10  ,7(1x,i5))
 
 
c  exits
c  -----
 
      return
900   call abt('ufbdmp - file is closed                     ')
901   call abt('ufbdmp - no message open                    ')
902   call abt('ufbdmp - i-node mismatch                    ')
      end
