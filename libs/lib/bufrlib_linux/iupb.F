      function iupb(mbay,nbyt,nbit)
 
      dimension mbay(*)
 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 
      mbit = (nbyt-1)*8
      call upb(iupb,nbit,mbay,mbit)
      return
      end
