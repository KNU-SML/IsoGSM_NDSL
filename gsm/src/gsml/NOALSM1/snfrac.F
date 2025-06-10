      subroutine snfrac (sneqv,snup,salp,snowh,sncovr)
c
      implicit none
c      
c ----------------------------------------------------------------------
c subroutine snfrac
c ----------------------------------------------------------------------
c calculate snow fraction (0 -> 1)
c sneqv   snow water equivalent (m)
c snup    threshold sneqv depth above which sncovr=1
c salp    tuning parameter
c sncovr  fractional snow cover
c ----------------------------------------------------------------------
      real sneqv, snup, salp, sncovr, rsnow, z0n, snowh
c ----------------------------------------------------------------------
c snup is veg-class dependent snowdepth threshhold (set in routine
c redprm) above which snocvr=1.
c ----------------------------------------------------------------------
          if (sneqv .lt. snup) then
            rsnow = sneqv/snup
            sncovr = 1. - ( exp(-salp*rsnow) - rsnow*exp(-salp))
          else
            sncovr = 1.0
          endif
c
          z0n=0.035 
c     formulation of dickinson et al. 1986
c        sncovr=snowh/(snowh + 5*z0n)
c     formulation of marshall et al. 1994
c        sncovr=sneqv/(sneqv + 2*z0n)
c ----------------------------------------------------------------------
c end subroutine snfrac
c ----------------------------------------------------------------------
      return
      end
