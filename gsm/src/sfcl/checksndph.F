      subroutine checksndph(snow,snwdph,ijdim)
c
      real snow(ijdim),snwdph(ijdim)
c
      do ij=1,ijdim
        if(snow(ij).gt.0.) then
          if(snwdph(ij).le.0.) then
            snwdph(ij)=snow(ij)*5.
          endif
        endif
      enddo
      return
      end
