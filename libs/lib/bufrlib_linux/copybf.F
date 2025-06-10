      subroutine copybf(lunin,lunot)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)

      character*8 sec0
      character*1 moct(24000)
      dimension   mbay(5000)
      equivalence (mbay(1),moct(1))

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      isec0 = 8/nbytw+1
      nmsg  = 0

c  check bufr file statuses
c  ------------------------

      call status(lunin,lun,il,im)
      if(il.ne.0) goto 900
      call status(lunot,lun,il,im)
      if(il.ne.0) goto 901

      rewind(lunin)
      rewind(lunot)

c  read and copy a bufr file on unit lunin to unit lunot
c  -----------------------------------------------------

1     read(lunin,end=2,err=902) sec0,(mbay(i),i=isec0,lmsg(sec0))
      write(lunot     ,err=903) sec0,(mbay(i),i=isec0,lmsg(sec0))
      goto 1

c  normal exit
c  -----------

2     close(lunin)
      close(lunot)
      return

c  error exits
c  -----------

900   call abt('copybf - input  file is currently open for bufr')
901   call abt('copybf - output file is currently open for bufr')
902   call abt('copybf - error reading file    ')
903   call abt('copybf - error writing file    ')
      end
