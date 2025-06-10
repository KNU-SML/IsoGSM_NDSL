      subroutine sfcfcnv(sfcftypin,sfcfcsin,ijdim,numsfcsin)
c
c  convert input surface sfcftypein file to model sfcftype 
c
c    sfcftypin: character ['osu1','osu2','noa1','vic1']
c    sfcfcsin : real input sfc file record
c    ijdim    : integer dimension of sfcfsin array
c    numsfcsin: integer number of input sfc file record
c
      implicit none
#include <paramodel.h>
c
      character*4 sfcftypin
      integer ijdim,numsfcsin
      real sfcfcsin(ijdim,numsfcsin)
c
      if(sfcftypin.eq.'osu1') then
c
#ifdef OSULSM1
         print *,'no surface file type conversion performed'
#endif
#ifdef OSULSM2
         call osu1toosu2(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef NOALSM1
         call osu1tonoa1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM1
         call osu1tovic1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM2
         call osu1tovic2(sfcfcsin,ijdim,numsfcsin)
#endif
c
      elseif(sfcftypin.eq.'osu2') then
c
#ifdef OSULSM1
         call osu2toosu1
#endif
#ifdef OSULSM2
         print *,'no surface file type conversion performed'
#endif
#ifdef NOALSM1
         call osu2tonoa1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM1
         call osu2tovic1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM2
         call osu2tovic2(sfcfcsin,ijdim,numsfcsin)
#endif
c
      elseif(sfcftypin.eq.'noa1') then
c
#ifdef OSULSM1
         call noa1toosu1
#endif
#ifdef OSULSM2
         call noa1toosu2
#endif
#ifdef NOALSM1
         print *,'no surface file type conversion performed'
#endif
#ifdef VICLSM1
         call noa1tovic1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM2
         call noa1tovic2(sfcfcsin,ijdim,numsfcsin)
#endif
c
      elseif(sfcftypin.eq.'vic1')then
c
#ifdef OSULSM1
         call vic1toosu1
#endif
#ifdef OSULSM2
         call vic1toosu2
#endif
#ifdef NOALSM1
         call vic1tonoa1
#endif
#ifdef VICLSM1
         print *,'no surface file type conversion performed'
#endif
#ifdef VICLSM2
         call vic1tovic2(sfcfcsin,ijdim,numsfcsin)
#endif
c
      elseif(sfcftypin.eq.'vic2')then
c
#ifdef OSULSM1
         call vic2toosu1
#endif
#ifdef OSULSM2
         call vic2toosu2
#endif
#ifdef NOALSM1
         call vic2tonoa1
#endif
#ifdef VICLSM1
         call vic2tovic1
#endif
#ifdef VICLSM2
         print *,'no surface file type conversion performed'
#endif
c
      elseif(sfcftypin.eq.'ncp1')then
c
#ifdef OSULSM1
         call ncp1tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu1
#endif
#ifdef OSULSM2
         call ncp1tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu2
#endif
#ifdef NOALSM1
         call ncp1tonoa1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM1
         call ncp1tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1tovic1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef VICLSM2
         call ncp1tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1tovic2(sfcfcsin,ijdim,numsfcsin)
#endif
      elseif(sfcftypin.eq.'ncp2')then
c
#ifdef OSULSM1
         call ncp2tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu1
#endif
#ifdef OSULSM2
         call ncp2tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu2
#endif
#ifdef NOALSM1
         call ncp2tonoa1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef OSULSM1
         call ncp2tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu1
#endif
#ifdef OSULSM2
         call ncp2tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu2
#endif
      elseif(sfcftypin.eq.'ncp3')then
c
#ifdef OSULSM1
         call ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu1
#endif
#ifdef OSULSM2
         call ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu2
#endif
#ifdef NOALSM1
         call ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
#endif
#ifdef OSULSM1
         call ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu1
#endif
#ifdef OSULSM2
         call ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
         call noa1toosu2
#endif
      else
         print *,'illegal intput surface file type'
         call abort
      endif
c
      return
      end
