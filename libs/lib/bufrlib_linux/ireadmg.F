      function ireadmg(lunit,subset,idate)
      character*8 subset
      call readmg(lunit,subset,idate,ireadmg)
      return

      entry ireadsb(lunit)
      call readsb(lunit,ireadsb)
      return
      end
