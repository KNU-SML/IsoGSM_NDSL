      subroutine readdx(lunit,lun,lundx)

      common /quiet/ iprt

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  get the bufr status of unit lundx
c  ---------------------------------

      call status(lundx,lud,ildx,imdx)

c  read a dx-table from the indicated source
c  -----------------------------------------

      if (lunit.eq.lundx) then
         if(iprt.ge.1) print100,lundx,lunit
         rewind lunit
         call rdbfdx(lunit,lun)
      elseif(ildx.ne.0) then
         if(iprt.ge.1) print101,lundx,lunit
         call cpbfdx(lud,lun)
      elseif(ildx.eq.0) then
         if(iprt.ge.1) print102,lundx,lunit
         rewind lundx
         call rdusdx(lundx,lun)
      else
         call abt('readdx - screwup')
      endif

100   format(' reading bufr dx-tables from ',i2,' to ',i2)
101   format(' copying bufr dx-tables from ',i2,' to ',i2)
102   format(' reading user dx-tables from ',i2,' to ',i2)

      return
      end
