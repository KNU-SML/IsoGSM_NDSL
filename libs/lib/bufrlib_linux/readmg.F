      subroutine readmg(lunit,subset,jdate,iret)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)

      character*8 sec0,subset
      character*1 tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      iret = 0

c  check the file status
c  ---------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900
      if(il.gt.0) goto 901
      call wtstat(lunit,lun,il, 1)

c  read a message into a message buffer
c  ------------------------------------

1     mbit = 0
      imsg = 8/nbytw+1
      read(lunit,err=902,end=100) sec0,(mbay(i,lun),i=imsg,lmsg(sec0))
c     call pkc(sec0,8,mbay(1,lun),mbit)
      call chrtrn(mbay(1,lun),sec0,8)

c  section 1
c  ---------

      iad1 = 8
      len1 = iupb(mbay(1,lun),iad1+ 1,24)
      len2 = iupb(mbay(1,lun),iad1+ 8, 1)
      mtyp = iupb(mbay(1,lun),iad1+ 9, 8)
      mear = iupb(mbay(1,lun),iad1+13, 8)
      mmon = iupb(mbay(1,lun),iad1+14, 8)
      mday = iupb(mbay(1,lun),iad1+15, 8)
      mour = iupb(mbay(1,lun),iad1+16, 8)

      jdate = mear*10**6 + mmon*10**4 + mday*10**2 + mour

      if(mtyp.eq.11) goto 1

c  section 2
c  ---------

      iad2 = iad1+len1
      len2 = iupb(mbay(1,lun),iad2+1,24) * len2

c  section 3
c  ---------

      iad3 = iad2+len2
      len3 = iupb(mbay(1,lun),iad3+1 ,24)
      jsub = iupb(mbay(1,lun),iad3+5 ,16)
      isub = iupb(mbay(1,lun),iad3+10,16)

      call numtab(lun,isub,subset,tab,itab)
      call nemtba(lun,subset,msgt,mstb,inod)
      if(msgt.ne.mtyp) goto 903

c  section 4
c  ---------

      iad4 = iad3+len3
      len4 = iupb(mbay(1,lun),iad4+1,24)

      mbyt(lun) = iad4+4

c  normal exit
c  -----------

      if(inod.ne.inode(lun)) call strcln
      idate(lun) = jdate
      inode(lun) = inod
      msub (lun) = jsub
      nsub (lun) = 0
      nmsg (lun) = nmsg(lun)+1
      iret = 0
      return

c  eof on attempted read
c  ---------------------

100   call closbf(lunit)
      inode(lun) = 0
      idate(lun) = 0
      subset = ' '
      jdate = 0
      iret = -1
      return

c  error exits
c  -----------

900   call abt('readmg - file is closed                   ')
901   call abt('readmg - file is open for output          ')
902   call abt('readmg - i/o error reading message        ')
903   call abt('readmg - msgtype mismatch for '//subset    )
      end
