      subroutine ipkm(cbay,nbyt,n)
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      character*8 cbay,cint
      equivalence(cint,int)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      if(nbyt.gt.nbytw) call abt('ipkm - nbyt>wrd len')
      int = irev(ishft(n,(nbytw-nbyt)*8))
      do i=1,nbyt
      cbay(i:i) = cint(i:i)
      enddo
      return
      end
