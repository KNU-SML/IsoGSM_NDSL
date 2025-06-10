      subroutine fld2msk(data,ijdim,cval,cond,out)
#include <paramodel.h>
c
      implicit none
      integer ijdim,ij
      real cval
      real data(ijdim),out(ijdim)
      character*2 cond
c
      if(cond.eq.'ge') then
        do ij=1,ijdim
          if(data(ij).ge.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      elseif(cond.eq.'gt') then
        do ij=1,ijdim
          if(data(ij).gt.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      elseif(cond.eq.'le') then
        do ij=1,ijdim
          if(data(ij).le.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      elseif(cond.eq.'lt') then
        do ij=1,ijdim
          if(data(ij).lt.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      elseif(cond.eq.'eq') then
        do ij=1,ijdim
          if(data(ij).eq.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      elseif(cond.eq.'ne') then
        do ij=1,ijdim
          if(data(ij).ne.cval) then
            out(ij)=1.
          else
            out(ij)=0.
          endif
        enddo
      else
        write(6,*) ' illegal cond in fld2msk.  cond=',cond
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
      return
      end
