      function ifbget(lunit)
 
      common /msgcwd/ nmsg(10),nsub(10),msub(10),inode(10),idate(10)
 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
 
c  make sure a file/message is open for input
c  ------------------------------------------
 
      call status(lunit,lun,il,im)
      if(il.ge.0) goto 900
      if(im.eq.0) goto 901
 
c  see if there is another subset in the message
c  ---------------------------------------------
 
      if(nsub(lun).lt.msub(lun)) then
         ifbget = 0
      else
         ifbget = -1
      endif
 
c  exit one way or another
c  -----------------------
 
      return
900   call abt('ifbget - file not open for input')
901   call abt('ifbget - no message open        ')
      end
