      subroutine writsa(lunxx,msgt,msgl)

      common /bufrmg/ msglen,msgtxt(5000)

      dimension msgt(*)

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      lunit = abs(lunxx)

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.lt.0) goto 901
      if(im.eq.0) goto 902

c  see if a memory message is waiting or forced
c  --------------------------------------------

      if(lunxx.lt.0) call closmg(lunit)

      if(msglen.gt.0) then
         msgl = msglen
         do n=1,msgl
         msgt(n) = msgtxt(n)
         enddo
         msglen = 0
      else
         msgl = 0
      endif

      if(lunxx.lt.0) return

c  pack up the subset and put it into the message
c  ----------------------------------------------

      call wrtree(lun)
      call msgupd(lunit,lun)

      return
900   call abt('writsa - file is closed                     ')
901   call abt('writsa - file is open for input             ')
902   call abt('writsa - no message open                    ')
      end
