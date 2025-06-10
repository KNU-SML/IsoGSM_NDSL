      subroutine pkb(nval,nbits,ibay,ibit)
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      dimension ibay(*)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      nwd  = ibit/nbitw + 1
      nbt  = mod(ibit,nbitw)
      ival = nval
      if(ishft(ival,-nbits).gt.0) ival = -1
      int = ishft(ival,nbitw-nbits)
      int = ishft(int,-nbt)
      msk = ishft(  -1,nbitw-nbits)
      msk = ishft(msk,-nbt)
      ibay(nwd) = irev(ior(iand(irev(ibay(nwd)),not(msk)),int))
      if(nbt+nbits.gt.nbitw) then
         int = ishft(ival,2*nbitw-(nbt+nbits))
         msk = ishft(  -1,2*nbitw-(nbt+nbits))
         ibay(nwd+1) = irev(ior(iand(irev(ibay(nwd+1)),not(msk)),int))
      endif
 
      ibit = ibit + nbits
 
      return
      end
