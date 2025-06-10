      subroutine ufbinx(lunit,imsg,isub,usr,i1,i2,iret,str)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

      character*(*) str
      character*10  cunit
      character*8   subset
      dimension     usr(i1,i2)
      logical       openit
      real*8        usr

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      write(cunit,'(a8,i2)') 'ON UNIT ',lunit

c  see if we need to open a file
c  -----------------------------

      call status(lunit,lun,il,im)
      openit = il.eq.0

      if(openit) then
         call openbf(lunit,'IN',lunit)
      else
         call rewnbf(lunit,0)
      endif

c  read subset isub from message imsg
c  ------------------------------------

      do i=1,imsg-1
      read(lunit,err=900,end=901) mstr
      enddo

      call readmg(lunit,subset,jdate,jret)
      if(jret.ne.0) goto 901

      do i=1,isub-1
      if(nsub(lun).gt.msub(lun)) goto 902
      ibit = mbyt(lun)*8
      call upb(nbyt,16,mbay(1,lun),ibit)
      mbyt(lun) = mbyt(lun) + nbyt
      nsub(lun) = nsub(lun) + 1
      enddo
      call readsb(lunit,jret)
      if(jret.ne.0) goto 902

      call ufbint(lunit,usr,i1,i2,iret,str)

c  leave the file as it was before
c  -------------------------------

      if(openit) then
         call closbf(lunit)
      else
         call rewnbf(lunit,1)
      endif

c  normal exit
c  -----------

      return
900   call abt('ufbinx - error reading messages '//cunit)
901   call abt('ufbinx - eof - reading messages '//cunit)
902   call abt('ufbinx - eom - reading subsets  '//cunit)
      end
