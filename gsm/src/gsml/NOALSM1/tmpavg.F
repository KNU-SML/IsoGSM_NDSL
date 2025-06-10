      subroutine tmpavg (tavg,tup,tm,tdn,zsoil,nsoil,k) 
c      
      implicit none
c      
c ----------------------------------------------------------------------
c subroutine tmpavg
c ----------------------------------------------------------------------
c calculate soil layer average temperature (tavg) in freezing/thawing
c layer using up, down, and middle layer temperatures (tup, tdn, tm),
c where tup is at top boundary of layer, tdn is at bottom boundary of
c layer.  tm is layer prognostic state temperature.
c ----------------------------------------------------------------------
      integer k
      integer nsoil
c
      real dz
      real dzh
      real t0
      real tavg
      real tdn
      real tm
      real tup
      real x0
      real xdn
      real xup
      real zsoil (nsoil)
c
      parameter(t0 = 2.7315e2)
c
      if (k .eq. 1) then
        dz = -zsoil(1)
      else
        dz = zsoil(k-1)-zsoil(k)
      endif
      dzh=dz*0.5
      if (tup .lt. t0) then
        if (tm .lt. t0) then
          if (tdn .lt. t0) then
c ----------------------------------------------------------------------
c tup, tm, tdn < t0
c ----------------------------------------------------------------------
            tavg = (tup + 2.0*tm + tdn)/ 4.0            
          else
c ----------------------------------------------------------------------
c tup & tm < t0,  tdn >= t0
c ----------------------------------------------------------------------
            x0 = (t0 - tm) * dzh / (tdn - tm)
            tavg = 0.5 * (tup*dzh+tm*(dzh+x0)+t0*(2.*dzh-x0)) / dz
          endif      
        else
          if (tdn .lt. t0) then
c ----------------------------------------------------------------------
c tup < t0, tm >= t0, tdn < t0
c ----------------------------------------------------------------------
            xup  = (t0-tup) * dzh / (tm-tup)
            xdn  = dzh - (t0-tm) * dzh / (tdn-tm)
            tavg = 0.5 * (tup*xup+t0*(2.*dz-xup-xdn)+tdn*xdn) / dz
          else
c ----------------------------------------------------------------------
c tup < t0, tm >= t0, tdn >= t0
c ----------------------------------------------------------------------
            xup  = (t0-tup) * dzh / (tm-tup)
            tavg = 0.5 * (tup*xup+t0*(2.*dz-xup)) / dz
          endif   
        endif
      else
        if (tm .lt. t0) then
          if (tdn .lt. t0) then
c ----------------------------------------------------------------------
c tup >= t0, tm < t0, tdn < t0
c ----------------------------------------------------------------------
            xup  = dzh - (t0-tup) * dzh / (tm-tup)
            tavg = 0.5 * (t0*(dz-xup)+tm*(dzh+xup)+tdn*dzh) / dz
          else
c ----------------------------------------------------------------------
c tup >= t0, tm < t0, tdn >= t0
c ----------------------------------------------------------------------
            xup  = dzh - (t0-tup) * dzh / (tm-tup)
            xdn  = (t0-tm) * dzh / (tdn-tm)
            tavg = 0.5 * (t0*(2.*dz-xup-xdn)+tm*(xup+xdn)) / dz
          endif   
        else
          if (tdn .lt. t0) then
c ----------------------------------------------------------------------
c tup >= t0, tm >= t0, tdn < t0
c ----------------------------------------------------------------------
            xdn  = dzh - (t0-tm) * dzh / (tdn-tm)
            tavg = (t0*(dz-xdn)+0.5*(t0+tdn)*xdn) / dz                 
          else
c ----------------------------------------------------------------------
c tup >= t0, tm >= t0, tdn >= t0
c ----------------------------------------------------------------------
            tavg = (tup + 2.0*tm + tdn) / 4.0
          endif
        endif
      endif
c ----------------------------------------------------------------------
c end subroutine tmpavg
c ----------------------------------------------------------------------
      return
      end
