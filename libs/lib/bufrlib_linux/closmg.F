      subroutine closmg(lunit)

      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.lt.0) goto 901
      if(im.ne.0) call msgwrt(lunit,mbay(1,lun),mbyt(lun))
      call wtstat(lunit,lun,il,0)

      return
900   call abt('closmg - file is closed            ')
901   call abt('closmg - file is open for input    ')
      end
