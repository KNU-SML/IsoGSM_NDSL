      subroutine openmb(lunit,subset,jdate)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)

      character*(*) subset
      logical       open

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.lt.0) goto 901

c  get some subset particulars
c  ---------------------------

      call nemtba(lun,subset,mtyp,mstb,inod)
      open = im.eq.0 .or. inod.ne.inode(lun) .or. jdate.ne.idate(lun)

c  maybe open a new or different type of message
c  ---------------------------------------------

      if(open) then
         call closmg(lunit)
         call wtstat(lunit,lun,il, 1)
         inode(lun) = inod
         idate(lun) = jdate
         call msgini(lun)
         call usrtpl(lun,1,1)
         call strcln
      endif

      return
900   call abt('openmb - file is closed            ')
901   call abt('openmb - file is open for input    ')
      end
