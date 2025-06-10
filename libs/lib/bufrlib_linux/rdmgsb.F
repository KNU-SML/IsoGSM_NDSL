      subroutine rdmgsb(lunit,imsg,isub)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

      character*8 subset

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  open the file and skip to message imsg
c  ----------------------------------------

      call openbf(lunit,'IN',lunit)
      call status(lunit,lun,il,im)

      do i=1,imsg-1
      read(lunit,err=900,end=901)
      enddo

      call readmg(lunit,subset,jdate,iret)
      if(iret.ne.0) goto 901

c  position at subset isub
c  -------------------------

      do i=1,isub-1
      if(nsub(lun).gt.msub(lun)) goto 902
      ibit = mbyt(lun)*8
      call upb(nbyt,16,mbay(1,lun),ibit)
      mbyt(lun) = mbyt(lun) + nbyt
      nsub(lun) = nsub(lun) + 1
      enddo

      call readsb(lunit,iret)
      if(iret.ne.0) goto 902

c  normal exit
c  -----------

      return
900   call abt('rdmgsb - error reading messages  ')
901   call abt('rdmgsb - eof - reading messages  ')
902   call abt('rdmgsb - eom - reading subsets   ')
      end
