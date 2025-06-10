      subroutine tbnd (tu,tb,zsoil,zbot,k,nsoil,tbnd1)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine tbnd
c ----------------------------------------------------------------------
c calculate temperature on the boundary of the layer by interpolation of
c the middle layer temperatures
c ----------------------------------------------------------------------
      integer nsoil
      integer k
c
      real tbnd1
      real t0
      real tu
      real tb
      real zb
      real zbot
      real zup
      real zsoil (nsoil)
c
      parameter(t0 = 273.15)
c ----------------------------------------------------------------------
c use surface temperature on the top of the first layer
c ----------------------------------------------------------------------
      if (k .eq. 1) then
        zup = 0.
      else
        zup = zsoil(k-1)
      endif
c ----------------------------------------------------------------------
c use depth of the constant bottom temperature when interpolate
c temperature into the last layer boundary
c ----------------------------------------------------------------------
      if (k .eq. nsoil) then
        zb = 2.*zbot-zsoil(k)
      else
        zb = zsoil(k+1)
      endif
c ----------------------------------------------------------------------
c linear interpolation between the average layer temperatures
c ----------------------------------------------------------------------
      tbnd1 = tu+(tb-tu)*(zup-zsoil(k))/(zup-zb)
c ----------------------------------------------------------------------
c end subroutine tbnd
c ----------------------------------------------------------------------
      return
      end
