      subroutine nfdr_nofuel(ifmpix,isc,ierc,ibi,ic)
c   
      implicit none
c
      integer ifmpix,isc,ierc,ibi,ic
c
c     If Ag land      
c
      if(ifmpix .eq. 13)then
      isc = 251
      ierc = 251
      ibi = 251
      ic = 251
      endif
c
c     If Barren
c
      if(ifmpix .eq. 23)then
      isc = 252
      ierc = 252
      ibi = 252
      ic = 252
      endif
c
c     If Marsh
c
      if(ifmpix .eq. 24)then
      isc = 253
      ierc = 253
      ibi = 253
      ic = 253
      endif
c
c     If outside US
c
      if(ifmpix .eq. 0)then
      isc = 254
      ierc = 254
      ibi = 254
      ic = 254
      endif
c
c     If water pixel
c
      if(ifmpix .eq. 22)then
      isc = 255
      ierc = 255
      ibi = 255
      ic = 255
      endif      
c
      return
      end     
