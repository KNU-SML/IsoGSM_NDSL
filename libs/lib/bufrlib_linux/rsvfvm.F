      subroutine rsvfvm(nem1,nem2)

      character*8 nem1,nem2

      do i=1,len(nem1)
      if(i.eq.1) then
         j = 1 
      else
         if(nem1(i:i).eq.'.') then
            nem1(i:i) = nem2(j:j)
            j = j+1    
         endif
      endif
      enddo

      return
      end
