      subroutine msgupd(lunit,lun)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /msgptr/ nby0,nby1,nby2,nby3,nby4,nby5
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

      character*8 cbay
      equivalence (cbay,jbay)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  pad the subset buffer
c  ---------------------

      call pad(ibay,ibit,ibyt,8)

c  see if the new subset fits
c  --------------------------

      if(mbyt(lun)+ibyt.gt.maxbyt) then
         call msgwrt(lunit,mbay(1,lun),mbyt(lun))
         call msgini(lun)
      endif

      if(ibyt.gt.maxbyt-mbyt(lun)) goto 900

c  set a byte count and transfer the subset buffer into the message
c  ----------------------------------------------------------------

      lbit = 0
      call pkb(ibyt,16,ibay,lbit)

      call mvb(ibay,1,mbay(1,lun),mbyt(lun)-3,ibyt)
c     lbit = (mbyt(lun)-4)*8
c     do i=1,ibyt/nbytw+1
c     call pkb(ibay(i),nbitw,mbay(1,lun),lbit)
c     enddo

c  update the subset and byte counters
c  --------------------------------------

      mbyt(lun)   = mbyt(lun)   + ibyt
      nsub(lun)   = nsub(lun)   + 1

      lbit = (nby0+nby1+nby2+4)*8
      call pkb(nsub(lun),16,mbay(1,lun),lbit)

      lbyt = nby0+nby1+nby2+nby3
      nbyt = iupb(mbay(1,lun),lbyt+1,24)
      lbit = lbyt*8
      call pkb(nbyt+ibyt,24,mbay(1,lun),lbit)

c  reset the user arrays
c  ---------------------

      call usrtpl(lun,1,1)

      return
900   call abt('msgupd - subset longer than any possible message')
      end
