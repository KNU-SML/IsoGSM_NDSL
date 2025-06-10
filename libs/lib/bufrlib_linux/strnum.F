      subroutine strnum(str,num)

      character*(*) str
      character*20  str2

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      num = 0
      k = 0

      call strsuc(str,str2,num)

      do i=1,num
      read(str(i:i),'(i1)',err=99) j
      if(j.eq.0 .and. str(i:i).ne.'0') goto 99
      k = k*10+j
      enddo

      num = k
      return

99    num = -1
      return
      end
