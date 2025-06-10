      subroutine nfdr_eqmc(t,h,equil)
c
      implicit none
c
      real t,h,equil
c
c     This function delivers equilibrium moisture content (percent) of
c     dead fuels given temperature (deg. f) and relative humidity
c     (percent), based upon regression of forest products lab tables.
c 
      if (h.gt.10.) go to 10
      equil=0.03229+0.281073*h-0.000578*h*t
      go to 30
   10 if (h.gt.50.) go to 20
      equil=2.22749+0.160107*h-0.014784*t
      go to 30
   20 equil=21.0606+0.005565*h**2-0.00035*h*t-0.483199*h
   30 continue
c
      return
      end
