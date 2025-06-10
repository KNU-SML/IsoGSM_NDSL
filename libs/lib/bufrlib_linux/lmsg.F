      function lmsg(sec0)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)

      character*8 sec0

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      imsg = 8/nbytw
      lmsg = iupm(sec0(5:7),24)/8
      if(lmsg.eq.0) return
      lmsg = (lmsg+1)*imsg
      return
      end
