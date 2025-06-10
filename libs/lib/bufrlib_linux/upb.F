      subroutine upb(nval,nbits,ibay,ibit)
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      dimension ibay(*)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      nwd = (ibit)/nbitw+1
      nbt = mod(ibit,nbitw)
      int = ishft(irev(ibay(nwd)),nbt)
      int = ishft(int,nbits-nbitw)
      lbt = nbt+nbits
      if(lbt.gt.nbitw) jnt = irev(ibay(nwd+1))
      if(lbt.gt.nbitw) int = ior(int,ishft(jnt,lbt-2*nbitw))
      ibit = ibit+nbits
      nval = int
      return
      end
