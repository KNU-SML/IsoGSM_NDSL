      subroutine acritn(len, pl, acr)
c
c  subroutine  acritn
c          01 february 1991
c
      dimension pl(len), acr(len)
      integer, allocatable :: iwk(:)
      dimension a(15), ac(15), ad(15)
      common /rasacr/ a, ac, ad, actop
c
      allocate (iwk(len))
c
      do i=1,len
         iwk(i) = pl(i) * 0.02 - 0.999999999
         if (iwk(i) .gt. 1) then
            if (iwk(i) .le. 15) then
               acr(i) = ac(iwk(i)) + pl(i) * ad(iwk(i))
            else
               acr(i) = a(15)
            endif
         else
            acr(i) = actop
         endif
      enddo
      deallocate (iwk)
c
      return
      end
