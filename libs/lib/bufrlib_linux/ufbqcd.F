      subroutine ufbqcd(lunit,nemo,qcd)

      character*(*) nemo
      character*6  fxy,adn30
      character*1  tab

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      call status(lunit,lun,il,im)
      if(il.eq.0) goto 900

      call nemtab(lun,nemo,idn,tab,iret)
      if(tab.ne.'D') goto 901

      fxy = adn30(idn,6)
      if(fxy(2:3).ne.'63') goto 902
      read(fxy(4:6),'(f3.0)',err=903) qcd

      return
900   call abt('ufbqcd - file is closed                       ')
901   call abt('ufbqcd - missing or invalid table d qc code   ')
902   call abt('ufbqcd - table d qc code descriptor not 363yyy')
903   call abt('ufbqcd - error reading yyy from qc code descrp')
      end
