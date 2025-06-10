      subroutine pkc(chr,nchr,ibay,ibit)
 
      common /charac/ iascii,iatoe(0:255),ietoa(0:255)
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      character*(*) chr
      character*1   cval(8)
      dimension     ibay(*)
      equivalence   (cval,ival)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      if(nchr.gt.len(chr)) call abt('pkc - chr < nchr')
      lb = iord(nbytw)
      ival = 0
      nbit = 8
 
      do i=1,nchr
      cval(lb) = chr(i:i)
      if(iascii.eq.0) call ipkm(cval(lb),1,ietoa(iupm(cval(lb),8)))
      nwd  = ibit/nbitw + 1
      nbt  = mod(ibit,nbitw)
      int = ishft(ival,nbitw-nbit)
      int = ishft(int,-nbt)
      msk = ishft(  -1,nbitw-nbit)
      msk = ishft(msk,-nbt)
      ibay(nwd) = irev(ior(iand(irev(ibay(nwd)),not(msk)),int))
      if(nbt+nbit.gt.nbitw) then
         int = ishft(ival,2*nbitw-(nbt+nbit))
         msk = ishft(  -1,2*nbitw-(nbt+nbit))
         ibay(nwd+1) = irev(ior(iand(irev(ibay(nwd+1)),not(msk)),int))
      endif
      ibit = ibit + nbit
      enddo
 
      return
      end
