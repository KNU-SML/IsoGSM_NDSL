      subroutine chrtrn(str,chr,n)

      character*(*) str
      character*1   chr(n)

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do i=1,n
      str(i:i) = chr(i)
      enddo
      return
      end
