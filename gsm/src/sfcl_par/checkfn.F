      subroutine checkfn(condir,bindir,fn,numgrbs,
     1                   fnmskg,fnorog,fnmask)
#include <paramodel.h>
c
c  checks sfc file names for consistency.
c
c  if condir and bin dir are given, add them in front
c  of fn, fntfc, fnmskg,fnorog and fnmask to
c  provide full directory structure.  
c
c  note that if filenames starts from '/', assume that
c  full directory is already given and condir, bindir are not added.
c
      implicit none
      integer numgrbs
      character*128 condir,bindir
      character*128 fn(numgrbs)
      character*128 fnmskg,fnorog,fnmask
c
      integer n,nbin,ncon
c
      n=1
      dowhile (condir(n:n).ne.' '.and.n.le.128)
        n=n+1
      enddo
      ncon=n-1
      n=1
      dowhile (bindir(n:n).ne.' '.and.n.le.128)
        n=n+1
      enddo
      nbin=n-1
c
      do n=1,numgrbs
        if(ncon.gt.0) then
          if(fn(n)(1:4).ne.'    '.and.fn(n)(1:1).ne.'/') then
             fn(n)=condir(1:ncon)//'/'//fn(n)
          endif
        endif
      enddo
      if(fnmskg(1:4).eq.'    ') then
        print *,'fnmskg empty'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      else
        if(ncon.gt.0.and.fnmskg(1:1).ne.'/') then
          fnmskg=condir(1:ncon)//'/'//fnmskg
        endif
      endif
      if(fnorog(1:4).eq.'    ') then
        print *,'fnorog empty'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      else
        if(nbin.gt.0.and.fnorog(1:1).ne.'/') then
          fnorog=bindir(1:nbin)//'/'//fnorog
        endif
      endif
      if(fnmask(1:4).eq.'    ') then
        print *,'fnmask empty'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      else
        if(nbin.gt.0.and.fnmask(1:1).ne.'/') then
          fnmask=bindir(1:nbin)//'/'//fnmask
        endif
      endif
c
      return
      end
