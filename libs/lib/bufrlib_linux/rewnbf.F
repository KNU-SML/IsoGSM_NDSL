      subroutine rewnbf(lunit,isr)

      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /bufrsr/ junn,jill,jimm,jbit,jbyt,jmsg,jsub,ksub,jnod,jdat,
     .                jsr(10),jbay(5000)

      character*20 mstr
      character*10 cunit

c     data jsr /10*0/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      write(cunit,'(a8,i2)') 'ON UNIT ',lunit

c  try to trap bad call problems
c  -----------------------------

      if(isr.eq.0) then
         call status(lunit,lun,il,im)
         if(jsr(lun).ne.0)  goto 900
         if(il.eq.0) goto 901
      elseif(isr.eq.1) then
         lun = junn
         if(jsr(junn).ne.1)  goto 902
      else
         goto 903
      endif

c  store file parameters and set for reading
c  -----------------------------------------

      if(isr.eq.0) then
         junn = lun
         jill = il
         jimm = im
         jbit = ibit
         jbyt = mbyt (lun)
         jmsg = nmsg (lun)
         jsub = nsub (lun)
         ksub = msub (lun)
         jnod = inode(lun)
         jdat = idate(lun)
         do i=1,jbyt
         jbay(i) = mbay(i,lun)
         enddo
         call wtstat(lunit,lun,-1,0)
      endif

c  rewind the file and position after the dictionary
c  -------------------------------------------------

      rewind lunit
1     read(lunit,err=904,end=2) mstr
      if(ichar(mstr(17:17)).eq.11) goto 1
2     backspace lunit

c  restore file parameters and position it to where it was saved
c  -------------------------------------------------------------

      if(isr.eq.1) then
         lun        = junn
         il         = jill
         im         = jimm
         ibit       = jbit
         mbyt (lun) = jbyt
         nmsg (lun) = jmsg
         nsub (lun) = jsub
         msub (lun) = ksub
         inode(lun) = jnod
         idate(lun) = jdat
         do i=1,jbyt
         mbay(i,lun) = jbay(i)
         enddo
         do imsg=1,jmsg
         read(lunit,err=904,end=905)
         enddo
         call wtstat(lunit,lun,il,im)
      endif

c  normal exit
c  -----------

      jsr(lun) = mod(jsr(lun)+1,2)

      return

c  error encountered
c  -----------------

900   call abt('rewnbf - attempt to save before restore '//cunit)
901   call abt('rewnbf - file is closed  '               //cunit)
902   call abt('rewnbf - attempt to restore before save '//cunit)
903   call abt('rewnbf - save/restore flag ne 0 or 1 '   //cunit)
904   call abt('rewnbf - file read error '               //cunit)
905   call abt('rewnbf - eof trying to repostion file '  //cunit)
      end
