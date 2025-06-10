      subroutine naopen(lu,cfn,bbuf,isize,iret)
c$$$  subprogram documentation block
c
c subprogram: naopen         byte-addressable open
c   prgmmr: kanamistu        date: 04-02-06
c
c abstract: open a file  and move entire grib record to an array
c           to be accessed by naread 
c
c program history log:
c   04-02-06  kanamitsu     standard f90 version
c
c usage:    call naopen(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character (*) file name to open
c   output arguments:
c     iret         integer return code (0 if successful)
c
c attributes:
c   language: fortran 90
c
c$$$
      implicit none
c
      integer lu,iret,isize
      character bbuf(isize)
      character cfn*(*)
c
      integer :: i
      integer :: ier
      integer :: iy
      integer :: status(15)
      integer :: fstat  ! service function to get file info. 
c
c reopen with the correct record length (=file size)
c
      open(lu,file=cfn,iostat=iret,access='direct',recl=isize,
     1     form='unformatted')
      read(lu,rec=1) bbuf
      close(lu)
c
c     print *,'bbuf=',(bbuf(i),i=1,100)
c
      return
      end
