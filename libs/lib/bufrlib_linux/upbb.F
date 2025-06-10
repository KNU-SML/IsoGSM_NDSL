      subroutine upbb(int,nbit,mbit,lun)
 
      common /bitbuf/ maxbyt,ibit,ibay(5000),mbyt(10),mbay(5000,10)
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      nwd = mbit/nbitw + 1
      nbt = mod(mbit,nbitw)
      lbt = nbt+nbit
      iba = irev(mbay(nwd,lun))
      int = ishft(ishft(iba,nbt),nbit-nbitw)
      if(lbt.gt.nbitw) then
         iba = irev(mbay(nwd+1,lun))
         int = ior(int,ishft(iba,lbt-2*nbitw))
      endif
      return
      end
