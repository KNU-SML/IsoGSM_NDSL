      subroutine parseq(str,tags,mtag,ntag)

      character*(*) str,tags(mtag)
      character*80  astr
      logical       word

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      astr = str
      lstr = len(str)
      ltag = len(tags(1))
      if(lstr.gt.80) goto 900
      ntag = 0
      nchr = 0
      word = .false.

      do 10 i=1,lstr

      if(.not.word .and. str(i:i).ne.' ') then
         ntag = ntag+1
         if(ntag.gt.mtag) goto 901
         tags(ntag) = ' '
      endif

      if(word .and. str(i:i).eq.' ') nchr = 0
      word = str(i:i).ne.' '

      if(word) then
         nchr = nchr+1
         if(nchr.gt.ltag) goto 902
         tags(ntag)(nchr:nchr) = str(i:i)
      endif

10    continue

      return
900   call abt('parseq - string too long  '//astr)
901   call abt('parseq - too many tags    '//astr)
902   call abt('parseq - tag is too long  '//astr)
      end
