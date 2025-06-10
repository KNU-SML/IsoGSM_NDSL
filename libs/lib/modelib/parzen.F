      subroutine parzen(ipoint, weight)
*
c compute parzen window weight
*
c input
c   the number of points
      integer ipoint
*
c output
c   weight
      real weight (*)
*
      do i = 1, ipoint
        weight (i) = 1.
     1   - abs((float(i) - 0.5*float(ipoint+1))/(0.5*float(ipoint+1)))
      end do
      return
      end
