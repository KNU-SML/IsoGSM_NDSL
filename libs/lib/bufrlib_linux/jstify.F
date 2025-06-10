      subroutine jstify

      character*(*) str
      character*1  sign

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      entry jstchr(str)

      lstr = len(str)

      if(str.eq.' ') goto 900
1     if(str(1:1).eq.' ') then
         str  = str(2:lstr)
         goto 1
      endif
      return

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      entry jstnum(str,sign,iret)

      iret = 0
      lstr = len(str)

      if(str.eq.' ') goto 900
2     if(str(1:1).eq.' ') then
         str  = str(2:lstr)
         goto 2
      endif
      if(str(1:1).eq.'+') then
         str  = str(2:lstr)
         sign = '+'
      elseif(str(1:1).eq.'-') then
         str  = str(2:lstr)
         sign = '-'
      else
         sign = '+'
      endif

      call strnum(str,num)
      if(num.lt.0) iret = -1
      return

900   call abt('jstify - blank string not allowed')
      end
