      subroutine datebf(lunit,iy,im,id,ih,idate)

      character*24  mstr

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      idate = -1

c  see if the file is already open to bufr interface (a no-no)
c  -----------------------------------------------------------

      call status(lunit,lun,jl,jm)
      if(jl.ne.0) call abt('datebf - file already open')

c  check for no bufr data or no data at all
c  ----------------------------------------

      rewind lunit
      read(lunit,end=100,err=100) mstr
      if(mstr(1:4).ne.'BUFR') goto 100

c  read to a data message and pick out the date
c  --------------------------------------------

1     read(lunit,end=100,err=100) mstr
      if(ichar(mstr(17:17)).eq.11) goto 1
c     iy = ichar(mstr(21:21))
      iy = mod(ichar(mstr(21:21)),100)
      im = ichar(mstr(22:22))
      id = ichar(mstr(23:23))
      ih = ichar(mstr(24:24))
      idate = iy*1e6+im*1e4+id*1e2+ih
100   return
      end
