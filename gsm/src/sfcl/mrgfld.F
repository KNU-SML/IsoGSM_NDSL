      subroutine mrgfld(sfcfcs,grbfld,
     1                  ijdim,kdim,
     2                  lsf,lflgrb,loismsk)
c
      implicit none
      integer ijdim,kdim
      real    sfcfcs(ijdim,kdim)
      real    grbfld(ijdim,kdim)
      integer loismsk(ijdim)
      logical lsf,lflgrb(5)
c
      logical  lall
      integer  i,ij,k
c
      if(.not.lsf) return
      lall=.true.
      do i=1,5
        if(.not.lflgrb(i)) lall=.false.
      enddo
c
      if(lall) then
        do k=1,kdim
          do ij=1,ijdim
             sfcfcs(ij,k)=grbfld(ij,k)
          enddo
        enddo
        return
      endif
c
      do i=1,5
        if(lflgrb(i)) then
          do k=1,kdim
            do ij=1,ijdim
              if(loismsk(ij).eq.i-1) then
                sfcfcs(ij,k)=grbfld(ij,k)
              endif
            enddo
          enddo
        endif
      enddo
c
      return
      end
