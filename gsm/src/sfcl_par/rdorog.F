      subroutine rdorog(lugb,fnorog,ijdim,slmask,orog)
c
c  read model orography
c
      implicit none
#include <paramodel.h>
c
      integer lugb,ijdim
      character*132 fnorog
      real slmask(ijdim)
      real orog(ijdim)
c
      integer, allocatable :: loismsk(:)
c
      integer ij,nch
c
c  max/min over ocean, no snow land, no-snow sea ice, snow covered land
c  and snow covered sea ice.  at this stage, slmask contain only 0 or 1.
c
      real vmaxmin(2,5)
      data vmaxmin/9000.,-5000.,
     1             9000.,-5000.,
     2                0.,    0.,
     3                0.,    0.,
     4                0.,    0.
     5             /
c
      if(fnorog(1:6).eq.'      ') then
         print *,'orography file name empty'
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
      call numchar(fnorog,nch)
      open(unit=lugb,file=fnorog(1:nch),status='old',form='unformatted',
     1     err=910)
      go to 911
  910   continue
        print *,'error in opening file ',fnorog(1:nch)
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
  911 continue
      write(6,*) ' file ',fnorog(1:nch),' opened. unit=',lugb
      read(lugb) orog
c
      allocate (loismsk(ijdim))
      do ij=1,ijdim
        loismsk(ij)=nint(slmask(ij))
      enddo
c
c  quality control by checking max/min over various surfaces
c
      call qcmxmn(orog,ijdim,loismsk,vmaxmin,'orog    ')
c
      deallocate (loismsk)
      close(lugb)
      return
      end
