      subroutine ufbqcp(lunit,qcp,nemo)

      character*(*) nemo
      character*6  fxy,adn30
      character*1  tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900

      idn = ifxy('363000')+ifix(qcp)
      call numtab(lun,idn,nemo,tab,iret)

      return
900   call abt('ufbqcp - file is closed                       ')
      end
