      subroutine readsb(lunit,iret)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

c-----------------------------------------------------------------------
cfpp$ expand(status)
c-----------------------------------------------------------------------

      iret = 0

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.gt.0) goto 901
      if(im.eq.0) then
         iret = -1
         return
      endif

c  see if there is another subset in the message
c  ---------------------------------------------

      if(nsub(lun).eq.msub(lun)) then
         iret = -1
         return
      endif

c  read the next subset and reset the pointers
c  -------------------------------------------

      ibit = mbyt(lun)*8
      call upb(nbyt,16,mbay(1,lun),ibit)
      call rdtree(lun)
      mbyt(lun) = mbyt(lun) + nbyt
      nsub(lun) = nsub(lun) + 1

      return
900   call abt('readsb - file is closed           ')
901   call abt('readsb - file is open for output  ')
902   call abt('readsb - no message open          ')
      end
