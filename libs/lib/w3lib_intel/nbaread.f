      subroutine nbaread(bbuf,isize,ib,nb,ka,a)
c$$$  subprogram documentation block
c
c subprogram: nbaread         byte-addressable read
c   prgmmr: kanamitsu         date: 04-02-06
c
c abstract: read a given number of bytes from an unblocked file,
c   skipping a given number of bytes.
c
c program history log:
c
c usage:    call baread(bbuf,isize,ib,nb,ka,a)
c   input arguments:
c     bbuf         grib record
c     isize        size of grib record
c     ib           integer number of bytes to skip
c     nb           integer number of bytes to read
c   output arguments:
c     ka           integer number of bytes actually read
c     a            character*1 (nb) buffer read
c
c attributes:
c   language: fortran 90
c   machine: workstations
c
c$$$
      character bbuf(*)
      character a(nb)
c
c
      if(ib.gt.isize) then
        ka=0
        return
      endif
      ka=min(nb,isize)
      do i=1,ka
        if(i+ib-1.gt.isize) return
        a(i)=bbuf(i+ib-1)
c       a(i)=bbuf(i+ib)
      enddo
c     print *,'nbaread: isize,ib,nb,ka=',isize,ib,nb,ka
c
      return
      end
