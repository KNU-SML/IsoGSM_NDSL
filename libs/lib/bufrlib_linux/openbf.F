      subroutine openbf(lunit,io,lundx)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /stbfr / iolun(10),iomsg(10)

      character*(*) io
      character*20  mstr
      character*4   bufr
      logical       skipdx,append

      data ifirst/0/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      if(ifirst.eq.0) then
         call wrdlen
         call bfrini
         ifirst = 1
      endif

c  see if a file can be opened
c  ---------------------------

      call status(lunit,lun,il,im)
      if(lun.eq.0) goto 900
      if(il .ne.0) goto 901

c  check for no bufr data or no data at all in an "in" file
c  --------------------------------------------------------

      if(io.eq.'IN' .and. lunit.eq.lundx) then
         rewind lunit
         read(lunit,end=100,err=902) mstr
         ibit = 0
         call upc(bufr,4,mstr,ibit)
         if(bufr.ne.'BUFR') goto 902
      endif

c  set initial open defaults
c  -------------------------

      rewind lunit
      nmsg (lun) = 0
      nsub (lun) = 0
      msub (lun) = 0
      inode(lun) = 0
      idate(lun) = 0
      skipdx = .false.
      append = .false.

c  decide how to setup the dictionary
c  ----------------------------------

      if(io.eq.'IN') then
         call wtstat(lunit,lun,-1,0)
         call readdx(lunit,lun,lundx)
      else if(io.eq.'OUT') then
         call wtstat(lunit,lun, 1,0)
         call writdx(lunit,lun,lundx)
      else if(io.eq.'APN') then
         call wtstat(lunit,lun, 1,0)
         call readdx(lunit,lun,lundx)
         call posapn(lunit)
      else
         goto 903
      endif

      return

c  file opened for input is empty - let readmg give the bad news
c  -------------------------------------------------------------

100   rewind lunit
      call wtstat(lunit,lun,-1,0)
      call dxinit(lun,0)
      return

c  error exits
c  -----------

900   call abt('openbf - too many files opened already       ')
901   call abt('openbf - file already open                   ')
902   call abt('openbf - input file has non-bufr data        ')
903   call abt('openbf - io must be one of "IN" "OUT" "APN"  ')
      end
