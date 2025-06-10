      function iupm(cbay,nbits)
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      character*8 cbay,cint
      equivalence(cint,int)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      if(nbits.gt.nbitw) call abt('iupm - nbits>wrd len')
      cint = cbay
      int  = irev(int)
      iupm = ishft(int,nbits-nbitw)
      return
      end
