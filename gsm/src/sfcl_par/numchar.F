      subroutine numchar(char,num)
c
      implicit none
      character*128 char
      integer num,n
c
      num=1
      dowhile(char(num:num).ne.' ')
        num=num+1
      enddo
      num=num-1
      do n=num+1,128
        char(n:n)=' '
      enddo
      return
      end
