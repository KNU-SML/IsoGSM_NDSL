      subroutine ngetgir(bbuf,isize,msk1,msk2,mnum,mbuf,cbuf,nlen,
     1                  nnum,iret)
c$$$  subprogram documentation block
c
c subprogram: getgir         reads a grib index file
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: read a grib file and return its index contents.
c   see subprogram ixgb for documentation of the index buffer structure.
c
c program history log:
c   95-10-31  iredell
c
c usage:    call getgir(bbuf,isize,msk1,msk2,mbuf,cbuf,nlen,nnum,iret)
c   input arguments:
c     bbuf         grib record
c     isize        size of grib record
c     msk1         integer number of bytes to search for first message
c     msk2         integer number of bytes to search for other messages
c     mnum         integer number of index records to skip (usually 0)
c     mbuf         integer length of cbuf in bytes
c   output arguments:
c     cbuf         character*1 (mbuf) buffer to receive index data
c     nlen         integer length of each index record in bytes
c     nnum         integer number of index records
c                  (=0 if no grib messages are found)
c     iret         integer return code
c                    0      all ok
c                    1      cbuf too small to hold index data
c
c subprograms called:
c   skgb           seek next grib message
c   ixgb           make index record
c
c remarks: subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      character bbuf(isize)
      character cbuf(mbuf)
      parameter(lindex=152)
c
c  search for first grib message
c
      iseek=0
      call nskgb(bbuf,isize,iseek,msk1,lskip,lgrib)
      do m=1,mnum
        if(lgrib.gt.0) then
          iseek=lskip+lgrib
          call nskgb(buf,isize,iseek,msk2,lskip,lgrib)
        endif
      enddo
c
c  make an index record for every grib record found
c
      nlen=lindex
      nnum=0
      iret=0
      dowhile(iret.eq.0.and.lgrib.gt.0)
        if(nlen*(nnum+1).le.mbuf) then
          nnum=nnum+1
          call nixgb(bbuf,isize,lskip,lgrib,nlen,nnum,cbuf)
          iseek=lskip+lgrib
          call nskgb(bbuf,isize,iseek,msk2,lskip,lgrib)
        else
          iret=1
        endif
      enddo
c
      return
      end
