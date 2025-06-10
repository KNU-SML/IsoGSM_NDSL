      subroutine chrtrna(str,chr,n)

      common /charac/ iascii,iatoe(0:255),ietoa(0:255)

      character*(*) str
      character*1   chr(n)

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do i=1,n
      str(i:i) = chr(i)
      if(iascii.eq.0) call ipkm(str(i:i),1,iatoe(iupm(str(i:i),8)))
      enddo
      return
      end
