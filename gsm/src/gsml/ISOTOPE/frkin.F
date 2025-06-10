#include <define.h>
      function frkin(kiso,wind)
c
c      surface kinetic fractionation by Merlivat and Jouzel 1979
c
      real frkin,wind
      integer kiso
c     
#ifdef NOFRAC
      frkin=1.0D0
#else
      if (kiso.eq.1) then       !! for 18O
         if (wind.ge.7.) then   
            frkin=1-(0.000285*wind+0.00082)
         else                   
            frkin=1-0.006
         endif
      elseif (kiso.eq.2) then   !! for D
         if (wind.ge.7.) then 
            frkin=1-(0.000285*wind+0.00082)*0.88
         else                 
            frkin=1-0.006*0.88
         endif
      endif
#endif
c
      return
      end

