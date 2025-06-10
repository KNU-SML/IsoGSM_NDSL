      function irev(n)
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      character*8 cint,dint
      equivalence(cint,int)
      equivalence(dint,jnt)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      if(nrev.eq.0) then
         irev = n
      else
         int = n
         do i=1,nbytw
         dint(i:i) = cint(iord(i):iord(i))
         enddo
         irev = jnt
      endif
 
      return
      end
