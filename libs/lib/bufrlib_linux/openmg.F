      subroutine openmg(lunit,subset,jdate)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)

      character*(*) subset

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.lt.0) goto 901
      if(im.ne.0) call closmg(lunit)
      call wtstat(lunit,lun,il, 1)

c  get some subset particulars
c  ---------------------------

      call nemtba(lun,subset,mtyp,mstb,inod)
      if(inode(lun).ne.inod) call strcln
      inode(lun) = inod
      idate(lun) = jdate

c  initialize the open message
c  ---------------------------

      call msgini(lun)
      call usrtpl(lun,1,1)

      return
900   call abt('openmg - file is closed            ')
901   call abt('openmg - file is open for input    ')
      end
