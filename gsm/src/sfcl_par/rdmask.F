      subroutine rdmask(lugb,fnmask,ijdim,slmask)
c
      implicit none
#include <paramodel.h>
c
c  land/sea mask
c
      integer lugb,ijdim,nch
      character*132 fnmask
      real slmask(ijdim)
c
      integer ij
c
      if(fnmask(1:6).eq.'      ') then
         print *,'land-sea mask file name empty'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
c
      close(lugb)
#ifdef ASSIGN
      call assign('assign -R')
#endif
      call numchar(fnmask,nch)
      open(unit=lugb,file=fnmask(1:nch),status='old',form='unformatted',
     1     err=920)
      go to 921
  920   continue
        print *,'error in opening file ',fnmask(1:nch)
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
  921 continue
      write(6,*) ' file ',fnmask(1:nch),' opened. unit=',lugb
      read(lugb) slmask
      do ij=1,ijdim
        slmask(ij)=nint(slmask(ij))
      enddo
c
      close(lugb)
c
      return
      end
