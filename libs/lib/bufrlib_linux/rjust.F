      function rjust(str)
      character*(*) str
      rjust = 0
      if(str.eq.' ') return
      lstr = len(str)
      do while(str(lstr:lstr).eq.' ')
         str = ' '//str(1:lstr-1)
      enddo
      return
      end
