      subroutine nskgb(bbuf,isize,iseek,mseek,lskip,lgrib)
c$$$  subprogram documentation block
c
c subprogram: skgb           search for next grib message
c   prgmmr: iredell          org: w/nmc23     date: 93-11-22
c
c abstract: this subprogram searches a file for the next grib 1 message.
c   a grib 1 message is identified by its indicator section, i.e.
c   an 8-byte sequence with 'grib' in bytes 1-4 and 1 in byte 8.
c   if found, the length of the message is decoded from bytes 5-7.
c   the search is done over a given section of the file.
c   the search is terminated if an eof or i/o error is encountered.
c
c program history log:
c   93-11-22  iredell
c   95-10-31  iredell   add call to nbaread 
c
c usage:    call nskgb(bbuf,isize,iseek,mseek,lskip,lgrib)
c   input arguments:
c     bbuf         grib record
c     isize        size of grib record
c     iseek        integer number of bytes to skip before search
c     mseek        integer maximum number of bytes to search
c   output arguments:
c     lskip        integer number of bytes to skip before message
c     lgrib        integer number of bytes in message (0 if not found)
c
c subprograms called:
c   nbaread       byte-addressable read
c   gbyte        get integer data from bytes
c
c attributes:
c   language: fortran
c
c$$$
      character bbuf(isize)
      parameter(lseek=128)
      character z(lseek)
c
      lgrib=0
      ks=iseek
      kn=min(lseek,mseek)
      kz=lseek
      dowhile(lgrib.eq.0.and.kn.ge.8.and.kz.eq.lseek)
        call nbaread(bbuf,isize,ks,kn,kz,z)
        kz=min(isize-ks+1,kz)
        km=kz-8+1
        k=0
        dowhile(lgrib.eq.0.and.k.lt.km)
          call gbyte(z,i4,(k+0)*8,4*8)
          call gbyte(z,i1,(k+7)*8,1*8)
          if(i4.eq.1196575042.and.i1.eq.1) then
            lskip=ks+k
            call gbyte(z,lgrib,(k+4)*8,3*8)
          endif
          k=k+1
        enddo
        ks=ks+km
        kn=min(lseek,iseek+mseek-ks)
      enddo
c
      return
      end
