      subroutine ufbcnt(lunit,kmsg,ksub)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check the file status - return the message and subset counters
c  --------------------------------------------------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) call abt('ufbcnt - file is closed')
      kmsg = nmsg(lun)
      ksub = nsub(lun)
      return
      end
