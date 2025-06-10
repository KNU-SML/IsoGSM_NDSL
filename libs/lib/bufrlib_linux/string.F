      subroutine string(str,lun,i1,io)

      parameter (jcons=52)

c-mk  common /usrstr/ jcon(52),jdummy(10)
      common /usrstr/ jcon(52)
      common /stcach/ mstr,nstr,lstr,luns(50),usrs(50),icon(52,50)

      character*(*) str
      character*80  usrs,ust

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      if(len(str).gt.80) goto 900
      ust  = str

c  see if string is in the cache
c  -----------------------------

      do i=1,nstr
      if(lun.eq.luns(i).and.str.eq.usrs(i)) then
         do j=1,jcons
         jcon(j) = icon(j,i)
         enddo
         if(jcon(1).gt.i1) goto 901
         return
      endif
      enddo

c  if not parse it and put it there
c  --------------------------------

c     print*,'str=',str
      call parusr(str,lun,i1,io)
      lstr = max(mod(lstr+1,mstr+1),1)
      nstr = min(nstr+1,mstr)
      luns(lstr) = lun
      usrs(lstr) = str
      do i=1,jcons
      icon(i,lstr) = jcon(i)
      enddo
      return
900   call abt('string - user string > 80 chars         :'//ust)
901   call abt('string - must be at least i1 store nodes:'//ust)
      end
