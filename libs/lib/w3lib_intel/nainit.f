      subroutine nainit(lu,cfn,isize,iret)
c$$$  subprogram documentation block
c
c subprogram: nainit         find grib file size
c   prgmmr: kanamistu        date: 04-02-06
c
c abstract: open a file  and get file size
c
c program history log:
c   04-02-06  kanamitsu     standard f90 version
c
c usage:    call naopen(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character (*) file name to open
c   output arguments:
c     isize        size of the file in byte
c     iret         integer return code (0 if successful)
c
c attributes:
c   language: fortran 90
c
c$$$
      implicit none
      character :: cfn*(*)
      integer :: lu, iret
c
      integer :: isize
      integer :: status(100)
      integer :: fstat  ! service function to get file info. 
c
      open(lu,file=cfn,iostat=iret,access='direct',recl=1,
     1     form='formatted')
c    1     form='unformatted')
      iret = fstat (lu, status )       ! get file info.
      if(iret.lt.0) then
         print *,'fstat error'
      endif
c
c allocate dimension
c
      print *,'status=',status
      write(6,*) 'file size=', status(8)
      isize=status(8)
c
      close( lu )
      return
      end
