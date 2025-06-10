      subroutine status(lunit,lun,il,im)

      common /stbfr/ iolun(10),iomsg(10)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      if(lunit.le.0 .or. lunit.gt.99) goto 900

c  clear the status indicators
c  ---------------------------

      lun = 0
      il  = 0
      im  = 0

c  see if the unit is defined
c  --------------------------

      do i=1,10
      if(abs(iolun(i)).eq.lunit) lun = i
      enddo

c  if not, check for file space - return lun=0 if no file space
c  ------------------------------------------------------------

      if(lun.eq.0) then
         do i=1,10
         if(iolun(i).eq.0) then
            lun = i
            return
         endif
         enddo
      endif

c  if file defined return statuses
c  -------------------------------

      il = sign(1,iolun(lun))
      im = iomsg(lun)

      return
900   call abt('status - illegal unit given')
      end
