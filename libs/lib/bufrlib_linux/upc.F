      subroutine upc(chr,nchr,ibay,ibit)
 
      common /charac/ iascii,iatoe(0:255),ietoa(0:255)
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
 
      character*(*) chr
      character*8   cval
      dimension     ibay(*)
      equivalence   (cval,ival)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      lb = iord(nbytw)
      do i=1,nchr
      call upb(ival,8,ibay,ibit)
      chr(i:i) = cval(lb:lb)
      if(iascii.eq.0) call ipkm(chr(i:i),1,iatoe(iupm(chr(i:i),8)))
      enddo
 
      return
      end
