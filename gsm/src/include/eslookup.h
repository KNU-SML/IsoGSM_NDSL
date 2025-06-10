c
c $Id: eslookup.h,v 1.3 2004/05/21 23:38:52 kana Exp $
c $Author: kana $
c
c
c Common block and statement functions for saturation vapor pressure
c look-up procedure, J. J. Hack, February 1990
c
      integer plenest  ! length of saturation vapor pressure table
      parameter (plenest=250)
c
c Table of saturation vapor pressure values es from tmin degrees
c to tmax+1 degrees k in one degree increments.  ttrice defines the
c transition region where es is a combination of ice & water values
c
      common/comes/estbl(plenest) ,tmin  ,tmax  ,ttrice ,pcf(6) ,
     $             epsqs          ,rgasv ,hlatf ,hlatv  ,cp     ,
     $             icephs
c
      real estbl      ! table values of saturation vapor pressure
      real tmin       ! min temperature (K) for table
      real tmax       ! max temperature (K) for table
      real ttrice     ! transition range from es over H2O to es over ice
      real pcf        ! polynomial coeffs -> es transition water to ice
      real epsqs      ! Ratio of h2o to dry air molecular weights 
      real rgasv      ! Gas constant for water vapor
      real hlatf      ! Latent heat of vaporization
      real hlatv      ! Latent heat of fusion
      real cp         ! specific heat of dry air
      logical icephs  ! false => saturation vapor press over water only
c
c Dummy variables for statement functions
c
      real td         ! dummy variable for function evaluation
      real tlim       ! intermediate variable for es look-up with estbl4
      real estblf     ! statement function es look-up
      real estbl4     ! statement function es look-up
c
c Statement functions used in saturation vapor pressure table lookup
c there are two ways to use these three statement functions.
c For compilers that do a simple in-line expansion:
c => ttemp = tlim(t)
c    es    = estbl4(ttemp)
c
c For compilers that provide real optimization:
c => es    = estblf(t)
c
      tlim(td) = max(min(td,tmax),tmin)
c
      estblf(td) =  (tmin + int(tlim(td)-tmin) - tlim(td) + 1.0)
     $            *estbl(int(tlim(td)-tmin)+1)
     $            -(tmin + int(tlim(td)-tmin) - tlim(td)      )
     $            *estbl(int(tlim(td)-tmin)+2)
c
      estbl4(td) =  (tmin+int(td-tmin)+1.0-td)*estbl(int(td-tmin)+1)
     $            + ( td-(tmin+int(td-tmin)) )*estbl(int(td-tmin)+2)
 
