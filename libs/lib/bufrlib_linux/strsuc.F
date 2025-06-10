      subroutine strsuc(str1,str2,lens)

      character*(*) str1,str2

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      lens = 0
      lstr = len(str1)

      do i=1,lstr
      if(str1(i:i).ne.' ') goto 2
      enddo
      return

2     do j=i,lstr
      if(str1(j:j).eq.' ') goto 3
      lens = lens+1
      str2(lens:lens) = str1(j:j)
      enddo
      return

3     do i=j,lstr
      if(str1(i:i).ne.' ') lens = -1
      enddo
      return

      end
