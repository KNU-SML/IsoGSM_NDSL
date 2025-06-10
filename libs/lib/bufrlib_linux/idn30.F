      function idn30(adn30,l30)

      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)

      character*(*) adn30

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      if(len(adn30).lt.l30) goto 900
      if(l30.eq.5) then
         read(adn30,'(i5)') idn30
         if(idn30.lt.0 .or. idn30.gt.65535) goto 901
      elseif(l30.eq.6) then
         idn30 = ifxy(adn30)
      else
         goto 902
      endif

      return
900   call abt('idn30 - function input string too short    ')
901   call abt('idn30 - idn out of range, not a descriptor ')
902   call abt('idn30 - character length l30 <> 5 or 6     ')
      end
