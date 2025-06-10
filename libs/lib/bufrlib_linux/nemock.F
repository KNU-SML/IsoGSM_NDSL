      function nemock(nemo)

      character*(*) nemo
      character*38  chrset

      data chrset /'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'/
      data nchr   /38/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  get the length of nemo
c  ----------------------

      lnemo = 0

      do i=len(nemo),1,-1
      if(nemo(i:i).ne.' ') then
         lnemo = i
         goto 1
      endif
      enddo

1     if(lnemo.lt.1 .or. lnemo.gt.8) then
         nemock = -1
         return
      endif

c  scan nemo for allowable characters
c  ----------------------------------

      do 10 i=1,lnemo
      do j=1,nchr
      if(nemo(i:i).eq.chrset(j:j)) goto 10
      enddo
      nemock = -1
      return
10    enddo

c  normal exit
c  -----------

      nemock = 0
      return
      end
