      subroutine basefnam(iy,im,id,ih,fh,fnam,nchar)
c
c  returns base field file name with verifying date as
c  a suffix (sig.yyyymmddhh).
c
      implicit none
      integer iy,im,id,ih,jy,jm,jd,jh,nchar
      real fh,rjday
      character*4 cy
      character*2 cm,cd,ch
      character*20 fnam
c
      nchar=15
c
      call vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
c
      if(jy.lt.1000) then
        print *,'year less than 1000'
        call abort
      endif
c
  100 format('0',i1)
  101 format(i2)
  102 format(i4)
      write(cy,102) jy
      if(jm.lt.10) then
        write(cm,100) jm
      else
        write(cm,101) jm
      endif
      if(jd.lt.10) then
        write(cd,100) jd
      else
        write(cd,101) jd
      endif
      if(jh.lt.10) then
        write(ch,100) jh
      else
        write(ch,101) jh
      endif
c
      fnam='base.'//cy//cm//cd//ch
c
      return
      end      
