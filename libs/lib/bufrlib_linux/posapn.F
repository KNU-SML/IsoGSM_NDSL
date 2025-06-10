      subroutine posapn(lunit)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)

      character*8 sec0
      dimension   mbay(5000)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  read and count messages in order to postion for append
c  ------------------------------------------------------

      imsg = 8/nbytw+1
      rewind lunit
      irec = 0

1     read(lunit,err=900,end=2) sec0,(mbay(i),i=imsg,lmsg(sec0))
      irec = irec+1
      goto 1

2     rewind lunit
      do j=1,irec
      read(lunit,err=900,end=901) sec0,(mbay(i),i=imsg,lmsg(sec0))
      enddo

      return
900   call abt('posapn - io err reading a message')
901   call abt('posapn - failure to read to eofle')
      end
