#define CP cp
#define G g
#define RD rd
#define RV rv
#define HFUS hfus
#define HVAP hvap

      module physconst
!---------------------------------------------------------
!---------------------------------------------------------
#include <paramodel.h>

        real, parameter :: rair   = rd_     ! Gas constant for dry air (J/K/kg) 
        real, parameter :: rh2o   = rv_     ! Gas constant for water vapor (J/K/kg)
        real, parameter :: cpair  = cp_     ! specific heat of dry air
        real, parameter :: gravit = g_      ! gravity
        real, parameter :: latice = hfus_   ! latent heat of fusion
        real, parameter :: latvap = hvap_   ! latent heat of vaporization
        real, parameter :: tmelt  = 273.16  ! Freezing point of water

        real, parameter :: epsilo = rair/rh2o ! ratio of h2o to dry air molecular weights
        real, parameter :: cappa  = rair/cpair ! R/Cp
        real, parameter :: zvir   = rh2o/rair - 1.  ! virtual temperature constant
        
      end module physconst
