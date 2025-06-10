      subroutine ufbrew(lunit,msg)

      character*20  mstr

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      msg = 0

c  make sure the file is closed to count messages
c  ----------------------------------------------

      call status(lunit,lun,il,im)
      if(il.ne.0) goto 900

c  rewind and count the data messages
c  ----------------------------------

      rewind lunit

1     read(lunit,end=100,err=901) mstr
      if(mstr(1:4).ne.'BUFR') goto 901
      if(ichar(mstr(17:17)).ne.11) msg = msg+1
      goto 1

100   close(lunit)
      return

c  error exits
c  -----------

900   call abt('ufbrew - file already open        ')
901   call abt('ufbrew - file has non-bufr data   ')
      end
