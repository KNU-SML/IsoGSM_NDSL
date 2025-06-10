#include <define.h>
      function freq(kiso,temp)
c
      real freq,temp
      integer kiso
c     
#ifdef NOFRAC
      freq=1.0D0
#else
      if (kiso.eq.1) then       !! for 18O
         if (temp.ge.273.) then !! evaporation 
            freq=
     $        exp( -0.002067D0 - 0.4156D0 / temp
     $        + 1137.D0 / (temp ** 2) )
         else                   !! sublimation
            freq=
     $        exp( -0.028224D0 + 11.839D0 / temp )
         endif
      elseif (kiso.eq.2) then   !! for D
         if (temp.ge.273.) then !! evaporation
            freq=
     $        exp( 0.052612D0 - 76.248D0 / temp
     $         + 24844.D0 / ( temp**2 ) )
         else
            freq=
     $        exp( -0.0945D0 + 16289.D0 / ( temp**2) )
         endif
      endif
#endif
c
      return
      end

