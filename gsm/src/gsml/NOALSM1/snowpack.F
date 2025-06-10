      subroutine snowpack (esd,dtsec,snowh,sndens,tsnow,tsoil)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine snowpack
c ----------------------------------------------------------------------
c calculate compaction of snowpack under conditions of increasing snow
c density, as obtained from an approximate solution of e. anderson's
c differential equation (3.29), noaa technical report nws 19, by victor
c koren, 03/25/95.
c ----------------------------------------------------------------------
c esd     water equivalent of snow (m)
c dtsec   time step (sec)
c snowh   snow depth (m)
c sndens  snow density (g/cm3=dimensionless fraction of h2o density)
c tsnow   snow surface temperature (k)
c tsoil   soil surface temperature (k)
c
c subroutine will return new values of snowh and sndens
c ----------------------------------------------------------------------
      integer ipol, j
c
      real bfac,c1,c2,sndens,dsx,dthr,dtsec,dw,snowhc,snowh,pexp,tavgc,
     &     tsnow,tsnowc,tsoil,tsoilc,esd,esdc,esdcx,g,kn
c
      parameter(c1 = 0.01, c2=21.0, g=9.81, kn=4000.0)
c
c ----------------------------------------------------------------------
c conversion into simulation units
c ----------------------------------------------------------------------
      snowhc = snowh*100.
      esdc = esd*100.
      dthr = dtsec/3600.
      tsnowc = tsnow-273.15
      tsoilc = tsoil-273.15
c ----------------------------------------------------------------------
c calculating of average temperature of snow pack
c ----------------------------------------------------------------------
      tavgc = 0.5*(tsnowc+tsoilc)                                    
c ----------------------------------------------------------------------
c calculating of snow depth and density as a result of compaction
c  sndens=ds0*(exp(bfac*esd)-1.)/(bfac*esd)
c  bfac=dthr*c1*exp(0.08*tavgc-c2*ds0)
c note: bfac*esd in sndens eqn above has to be carefully treated
c numerically below:
c   c1 is the fractional increase in density (1/(cm*hr)) 
c   c2 is a constant (cm3/g) kojima estimated as 21 cms/g
c ----------------------------------------------------------------------
      if (esdc .gt. 1.e-2) then
        esdcx = esdc
      else
        esdcx = 1.e-2
      endif
      bfac = dthr*c1*exp(0.08*tavgc-c2*sndens)
c     dsx = sndens*((dexp(bfac*esdc)-1.)/(bfac*esdc))
c ----------------------------------------------------------------------
c the function of the form (e**x-1)/x imbedded in above expression
c for dsx was causing numerical difficulties when the denominator "x"
c (i.e. bfac*esdc) became zero or approached zero (despite the fact that
c the analytical function (e**x-1)/x has a well defined limit as 
c "x" approaches zero), hence below we replace the (e**x-1)/x 
c expression with an equivalent, numerically well-behaved 
c polynomial expansion.
c
c number of terms of polynomial expansion, and hence its accuracy, 
c is governed by iteration limit "ipol".
c      ipol greater than 9 only makes a difference on double
c            precision (relative errors given in percent %).
c       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
c       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
c       ipol=7, for rel.error <~ 1.8 e-4 % ...
c ----------------------------------------------------------------------
      ipol = 4
      pexp = 0.
      do j = ipol,1,-1
c        pexp = (1. + pexp)*bfac*esdc/real(j+1) 
        pexp = (1. + pexp)*bfac*esdcx/real(j+1) 
      end do
      pexp = pexp + 1.
      dsx = sndens*(pexp)
c ----------------------------------------------------------------------
c above line ends polynomial substitution
c ----------------------------------------------------------------------
c     end of korean formulation
c
c     base formulation (cogley et al., 1990)
c     convert density from g/cm3 to kg/m3
c       dsm=sndens*1000.0
c 
c       dsx=dsm+dtsec*0.5*dsm*g*esd/
c    &      (1e7*exp(-0.02*dsm+kn/(tavgc+273.16)-14.643))
c 
c     convert density from kg/m3 to g/cm3
c       dsx=dsx/1000.0
c
c     end of cogley et al. formulation 
c
c ----------------------------------------------------------------------
c set upper/lower limit on snow density
c ----------------------------------------------------------------------
      if (dsx .gt. 0.40) dsx = 0.40
      if (dsx .lt. 0.05) dsx = 0.05
      sndens = dsx
c ----------------------------------------------------------------------
c update of snow depth and density depending on liquid water during
c snowmelt.  assumed that 13% of liquid water can be stored in snow per
c day during snowmelt till snow density 0.40.
c ----------------------------------------------------------------------
      if (tsnowc .ge. 0.) then
        dw = 0.13*dthr/24.
        sndens = sndens*(1.-dw)+dw
        if (sndens .gt. 0.40) sndens = 0.40
      endif
c ----------------------------------------------------------------------
c calculate snow depth (cm) from snow water equivalent and snow density.
c change snow depth units to meters
c ----------------------------------------------------------------------
      snowhc = esdc/sndens
      snowh = snowhc*0.01
c ----------------------------------------------------------------------
c end subroutine snowpack
c ----------------------------------------------------------------------
      return
      end
