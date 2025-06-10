      subroutine writsb(lunit)

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.lt.0) goto 901
      if(im.eq.0) goto 902

c  pack up the subset and put it into the message
c  ----------------------------------------------

      call wrtree(lun)
      call msgupd(lunit,lun)

      return
900   call abt('writsb - file is closed                     ')
901   call abt('writsb - file is open for input             ')
902   call abt('writsb - no message open                    ')
      end
