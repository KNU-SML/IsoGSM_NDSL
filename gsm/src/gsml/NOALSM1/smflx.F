      subroutine smflx (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,
     &                  sh2o,slope,kdt,frzfact,
     &                  smcmax,bexp,pc,smcwlt,dksat,dwsat,
     &                  smcref,shdfac,cmcmax,
     &                  smcdry,cfactr,runoff1,runoff2,runoff3,
     &                  edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,
     &                  drip)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine smflx
c ----------------------------------------------------------------------
c calculate soil moisture flux.  the soil moisture content (smc - a per
c unit volume measurement) is a dependent variable that is updated with
c prognostic eqns. the canopy moisture content (cmc) is also updated.
c frozen ground version:  new states added: sh2o, and frozen ground
c correction factor, frzfact and parameter slope.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer i
      integer k
      integer nsoil
      integer nroot
c
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
c
      real bexp
      real beta
      real cfactr
      real cmc
      real cmc2ms
      real cmcmax
      real devap
      real dksat
      real drip
      real dt
      real dummy
      real dwsat
      real ec
      real edir
      real et(nsoil)
      real eta1
      real etp1
      real ett
      real excess
      real frzfact
      real fxexp
      real kdt
      real pc
      real pcpdrp
      real prcp1
      real q2
      real rhsct
      real rhstt(nsold)
      real rtdis(nsoil)
      real runoff1
      real runoff2
      real runoff3
      real sfctmp
      real shdfac
      real smc(nsoil)
      real sh2o(nsoil)
      real sice(nsold)
      real sh2oa(nsold)
      real sh2ofg(nsold)
      real slope
      real smcdry
      real smcmax
      real smcref
      real smcwlt
      real trhsct
      real zsoil(nsoil)
c ----------------------------------------------------------------------
c temperature criteria for snowfall tfreez should have same value as in
c sflx.f (main subroutine)
c ----------------------------------------------------------------------
      real tfreez
      parameter(tfreez = 273.15)
c ----------------------------------------------------------------------
c executable code begins here if the potential evapotranspiration is
c greater than zero.
c ----------------------------------------------------------------------
      dummy = 0.
      edir = 0.
      ec = 0.
      ett = 0.
      do k = 1,nsoil
        et(k) = 0.
      end do
c
      if (etp1 .gt. 0.0) then
c ----------------------------------------------------------------------
c retrieve direct evaporation from soil surface.  call this function
c only if veg cover not complete.
c frozen ground version:  sh2o states replace smc states.
c ----------------------------------------------------------------------
        if (shdfac .lt. 1.) then
c** clu_rev: change sh2o(1) to smc(1)
!clu      edir = devap(etp1,sh2o(1),zsoil(1),shdfac,smcmax,
          edir = devap(etp1,smc(1),zsoil(1),shdfac,smcmax,
     &                 bexp,dksat,dwsat,smcdry,smcref,smcwlt,fxexp)
        endif
c ----------------------------------------------------------------------
c initialize plant total transpiration, retrieve plant transpiration,
c and accumulate it for all soil layers.
c ----------------------------------------------------------------------
        if (shdfac.gt.0.0) then
c
          call transp (et,nsoil,etp1,sh2o,cmc,zsoil,shdfac,smcwlt,
     &                 cmcmax,pc,cfactr,smcref,sfctmp,q2,nroot,rtdis)
c
          do k = 1,nsoil
            ett = ett + et ( k )
          end do
c ----------------------------------------------------------------------
c calculate canopy evaporation.
c if statements to avoid tangent linear problems near cmc=0.0.
c ----------------------------------------------------------------------
          if (cmc .gt. 0.0) then
            ec = shdfac * ( ( cmc / cmcmax ) ** cfactr ) * etp1
          else
            ec = 0.0
          endif
c ----------------------------------------------------------------------
c ec should be limited by the total amount of available water on the
c canopy.  -f.chen, 18-oct-1994
c ----------------------------------------------------------------------
          cmc2ms = cmc / dt
          ec = min ( cmc2ms, ec )
        endif
      endif
c ----------------------------------------------------------------------
c total up evap and transp types to obtain actual evapotransp
c ----------------------------------------------------------------------
      eta1 = edir + ett + ec
c ----------------------------------------------------------------------
c compute the right hand side of the canopy eqn term ( rhsct )
c ----------------------------------------------------------------------
      rhsct = shdfac * prcp1 - ec
c ----------------------------------------------------------------------
c convert rhsct (a rate) to trhsct (an amount) and add it to existing
c cmc.  if resulting amt exceeds max capacity, it becomes drip and will
c fall to the grnd.
c ----------------------------------------------------------------------
      drip = 0.
      trhsct = dt * rhsct
      excess = cmc + trhsct
      if (excess .gt. cmcmax) drip = excess - cmcmax
c ----------------------------------------------------------------------
c pcpdrp is the combined prcp1 and drip (from cmc) that goes into the
c soil
c ----------------------------------------------------------------------
      pcpdrp = (1. - shdfac) * prcp1 + drip / dt
c ----------------------------------------------------------------------
c store ice content at each soil layer before calling srt & sstep
c ----------------------------------------------------------------------
      do i = 1,nsoil
        sice(i) = smc(i) - sh2o(i)
      end do
c ----------------------------------------------------------------------
c call subroutines srt and sstep to solve the soil moisture
c tendency equations. 
c
c if the infiltrating precip rate is nontrivial,
c   (we consider nontrivial to be a precip total over the time step 
c    exceeding one one-thousandth of the water holding capacity of 
c    the first soil layer)
c then call the srt/sstep subroutine pair twice in the manner of 
c   time scheme "f" (implicit state, averaged coefficient)
c   of section 2 of kalnay and kanamitsu (1988, mwr, vol 116, 
c   pages 1945-1958)to minimize 2-delta-t oscillations in the 
c   soil moisture value of the top soil layer that can arise because
c   of the extreme nonlinear dependence of the soil hydraulic 
c   diffusivity coefficient and the hydraulic conductivity on the
c   soil moisture state
c otherwise call the srt/sstep subroutine pair once in the manner of
c   time scheme "d" (implicit state, explicit coefficient) 
c   of section 2 of kalnay and kanamitsu
c pcpdrp is units of kg/m**2/s or mm/s, zsoil is negative depth in m 
c ----------------------------------------------------------------------
c      if ( pcpdrp .gt. 0.0 ) then
      if ( (pcpdrp*dt) .gt. (0.001*1000.0*(-zsoil(1))*smcmax) ) then
c ----------------------------------------------------------------------
c frozen ground version:
c smc states replaced by sh2o states in srt subr.  sh2o & sice states
c included in sstep subr.  frozen ground correction factor, frzfact
c added.  all water balance calculations using unfrozen water
c ----------------------------------------------------------------------
        call srt (rhstt,edir,et,sh2o,sh2o,nsoil,pcpdrp,zsoil,
     &            dwsat,dksat,smcmax,bexp,runoff1, 
     &            runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
        call sstep (sh2ofg,sh2o,dummy,rhstt,rhsct,dt,nsoil,smcmax,
     &              cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
        do k = 1,nsoil
          sh2oa(k) = (sh2o(k) + sh2ofg(k)) * 0.5
        end do
        call srt (rhstt,edir,et,sh2o,sh2oa,nsoil,pcpdrp,zsoil,
     &            dwsat,dksat,smcmax,bexp,runoff1,
     &            runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
        call sstep (sh2o,sh2o,cmc,rhstt,rhsct,dt,nsoil,smcmax,
     &              cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
      else
        call srt (rhstt,edir,et,sh2o,sh2o,nsoil,pcpdrp,zsoil,
     &            dwsat,dksat,smcmax,bexp,runoff1,
     &            runoff2,dt,smcwlt,slope,kdt,frzfact,sice,ai,bi,ci)
        call sstep (sh2o,sh2o,cmc,rhstt,rhsct,dt,nsoil,smcmax,
     &              cmcmax,runoff3,zsoil,smc,sice,ai,bi,ci)
      endif
c      runof = runoff
c ----------------------------------------------------------------------
c end subroutine smflx
c ----------------------------------------------------------------------
      return
      end
