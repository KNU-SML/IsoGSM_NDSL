      subroutine snowz0 (sncovr,z0)
c
      implicit none
c      
c ----------------------------------------------------------------------
c subroutine snowz0
c ----------------------------------------------------------------------
c calculate total roughness length over snow
c sncovr  fractional snow cover
c z0      roughness length (m)
c z0s     snow roughness length:=0.001 (m)
c ----------------------------------------------------------------------
      real sncovr, z0, z0s
c      parameter (z0s=0.001)
c current noah lsm condition - mbek, 09-oct-2001
      z0s = z0
c
      z0 = (1-sncovr)*z0 + sncovr*z0s
c ----------------------------------------------------------------------
c end subroutine snowz0
c ----------------------------------------------------------------------
      return
      end
