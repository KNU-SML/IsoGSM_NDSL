      subroutine wtstat(lunit,lun,il,im)

      common /stbfr/ iolun(10),iomsg(10)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c  check on the arguments
c  ----------------------

      if(lunit.le.0)            goto 900
      if(lun  .le.0)            goto 901
      if(il.lt.-1 .or. il.gt.1) goto 902
      if(im.lt. 0 .or. il.gt.1) goto 903

c  check on lunit-lun combination
c  ------------------------------

      if(abs(iolun(lun)).ne.lunit) then
         if(iolun(lun).ne.0) goto 905
      endif

c  reset the file statuses
c  -----------------------

      if(il.ne.0) then
         iolun(lun) = sign(lunit,il)
         iomsg(lun) = im
      else
         iolun(lun) = 0
         iomsg(lun) = 0
      endif

      return
900   call abt('wtstat - bad lunit                               ')
901   call abt('wtstat - bad lun                                 ')
902   call abt('wtstat - bad il                                  ')
903   call abt('wtstat - bad im                                  ')
905   call abt('wtstat - attempt to redefine exisiting file unit ')
      end
