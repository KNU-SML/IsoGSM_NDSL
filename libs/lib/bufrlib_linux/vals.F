      function vals(str)

      character*(*) str
      character*99  astr
      character*8   fmt

      data noinline /0/

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      call strsuc(str,astr,lens)
      write(fmt,'(''(f'',i2,''.0)'')') lens
      read(str,fmt,err=900) val
      vals = val
      return
900   vals = 10e10
      return
      end
