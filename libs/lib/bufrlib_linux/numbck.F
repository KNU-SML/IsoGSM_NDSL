      function numbck(numb)

      character*10 chrset
      character*6  numb
      character*1  fc

      data chrset /'0123456789'/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      numbck = 0
      lnumb  = 0
      fc     = numb(1:1)

c  check the first character of numb
c  ---------------------------------

      if(.not.(fc.eq.'A' .or. fc.eq.'0' .or. fc.eq.'3')) goto 900

c  check the rest of numb
c  ----------------------

      do 10 i=2,6
      do j=1,10
      if(numb(i:i).eq.chrset(j:j)) goto 10
      enddo
      goto 900
10    enddo

c  check for a valid descriptor
c  ----------------------------

      read(numb,'(1x,i2,i3)',err=900) ix,iy
      if(ix.lt.0 .or. ix.gt. 63) goto 900
      if(iy.lt.0 .or. iy.gt.255) goto 900

c  normal exit
c  -----------

      numbck = 0
      return

c  error exit
c  ----------

900   numbck = -1
      return
      end
