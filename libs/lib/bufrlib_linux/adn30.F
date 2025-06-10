      function adn30(idn,l30)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)

      character*(*) adn30

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      if(len(adn30).lt.l30         ) goto 900
      if(idn.lt.0 .or. idn.gt.65535) goto 901
      if(l30.eq.5) then
         write(adn30,'(i5)') idn
      elseif(l30.eq.6) then
         idf = ishft(idn,-14)
         idx = ishft(ishft(idn,nbitw-14),-(nbitw-6))
         idy = ishft(ishft(idn,nbitw- 8),-(nbitw-8))
         write(adn30,'(i1,i2,i3)') idf,idx,idy
      else
         goto 902
      endif

      do i=1,l30
      if(adn30(i:i).eq.' ') adn30(i:i) = '0'
      enddo

      return
900   call abt('adn30 - function return string too short')
901   call abt('adn30 - idn out of range                ')
902   call abt('adn30 - character length l30 <> 5 or 6')
      end
