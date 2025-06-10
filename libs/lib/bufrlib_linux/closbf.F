      subroutine closbf(lunit)
      call status(lunit,lun,il,im)
      if(il.gt.0 .and. im.ne.0) call closmg(lunit)
      call wtstat(lunit,lun,0,0)
      close(lunit)
      return
      end
