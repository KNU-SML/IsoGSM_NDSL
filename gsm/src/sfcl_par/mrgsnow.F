#include <paramodel.h>
      subroutine mrgsnow(grbsnow,sfcsnow,sfctemp,
#ifdef NOALSM1
     1             sfcsnod,
#endif
     2             ijdim)
c
c   if snow cover available but no snow depth available
c     1. set fcst snow depth to zero where snow cover mask is 0
c     2. use fcst snow depth where snow cover mask is 1
c     3. compute snow depth where fcst snow cover mask is 0
c
      implicit none  
      real snwmax,snwmin
      parameter(snwmin=25.0,snwmax=100.)
c
      integer ijdim
      real grbsnow(ijdim),sfcsnow(ijdim),sfctemp(ijdim)
#ifdef NOALSM1
      real sfcsnod(ijdim)
#endif
      integer ij
c
      do ij=1,ijdim
c 1.
        if(grbsnow(ij).eq.0.) then
          sfcsnow(ij)=0.
#ifdef NOALSM1
          sfcsnod(ij)=0.
#endif
        endif
c 2.
c       if(grbsnow(ij).eq.1..and.sfcsnow(ij).gt.0.) then
c         sfcsnow(ij)=sfcsnow(ij)
c       endif
c 3.
        if(grbsnow(ij).eq.-999..and.sfcsnow(ij).eq.0.) then
          if(sfctemp(ij).lt.243.0) then
            sfcsnow(ij) = snwmax
          elseif(sfctemp(ij).lt.273.0) then
            sfcsnow(ij)=snwmin+
     1           (snwmax-snwmin)*(273.0-sfctemp(ij))/30.0
          else
            sfcsnow(ij)=snwmin
          endif
#ifdef NOALSM1
          sfcsnod(ij)=sfcsnow(ij)*5.
#endif
        endif
      enddo
c
      return
      end
