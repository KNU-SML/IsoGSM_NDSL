      subroutine snopac (etp,eta,prcp,prcp1,snowng,smc,smcmax,smcwlt,
     &                   smcref,smcdry,cmc,cmcmax,nsoil,dt,
     &                   sbeta,df1,
     &                   q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,stc,epsca,
     &                   sfcprs,bexp,pc,rch,rr,cfactr,sncovr,esd,sndens,
     &                   snowh,sh2o,slope,kdt,frzfact,psisat,snup,
     &                   zsoil,dwsat,dksat,tbot,zbot,shdfac,runoff1,
     &                   runoff2,runoff3,edir,ec,et,ett,nroot,snomlt,
     &                   ice,rtdis,quartz,fxexp,csoil,
     &                   beta,drip,dew,flx1,flx2,flx3)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine snopac
c ----------------------------------------------------------------------
c calculate soil moisture and heat flux values & update soil moisture
c content and soil heat content values for the case when a snow pack is
c present.
c ----------------------------------------------------------------------
      integer ice
      integer nroot
      integer nsoil
c
      logical snowng
c
      real bexp
      real beta
      real cfactr
      real cmc
      real cmcmax
      real cp
      real cph2o
      real cpice
      real csoil
      real denom
      real dew
      real df1
      real dksat
      real drip
      real dsoil
      real dtot
      real dt
      real dwsat
      real ec
      real edir
      real epsca
      real esd
      real esdmin
      real expsno
      real expsoi
      real eta
      real eta1
      real etp
      real etp1
      real etp2
      real et(nsoil)
      real ett
      real ex
      real expfac
      real fdown
      real fxexp
      real flx1
      real flx2
      real flx3
      real f1
      real kdt
      real lsubf
      real lsubc
      real lsubs
      real pc
      real prcp
      real prcp1
      real q2
      real rch
      real rr
      real rtdis(nsoil)
      real ssoil
      real sbeta
      real ssoil1
      real sfctmp
      real shdfac
      real sigma
      real smc(nsoil)
      real sh2o(nsoil)
      real smcdry
      real smcmax
      real smcref
      real smcwlt
      real snomlt
      real snowh
      real stc(nsoil)
      real t1
      real t11
      real t12
      real t12a
      real t12b
      real t24
      real tbot
      real zbot
      real th2
      real yy
      real zsoil(nsoil)
      real zz1
      real tfreez
      real salp
      real sfcprs
      real slope
      real frzfact
      real psisat
      real snup
      real runoff1
      real runoff2
      real runoff3
      real quartz
      real sndens
      real sncond
      real rsnow
      real sncovr
      real qsat
      real etp3
      real seh
      real t14
      real csnow
c
      parameter(cp = 1004.5)
      parameter(cph2o = 4.218e+3)
      parameter(cpice = 2.106e+3)
      parameter(esdmin = 1.e-6)
      parameter(lsubf = 3.335e+5)
      parameter(lsubc = 2.501000e+6)
      parameter(lsubs = 2.83e+6)
      parameter(sigma = 5.67e-8)
      parameter(tfreez = 273.15)
c ----------------------------------------------------------------------
c executable code begins here:
c convert potential evap (etp) from kg m-2 s-1 to m s-1 and then to an
c amount (m) given timestep (dt) and call it an effective snowpack
c reduction amount, etp2 (m).  this is the amount the snowpack would be
c reduced due to evaporation from the snow sfc during the timestep.
c evaporation will proceed at the potential rate unless the snow depth
c is less than the expected snowpack reduction.
c if seaice (ice=1), beta remains=1.
c ----------------------------------------------------------------------
      prcp1 = prcp1*0.001
c
      etp2 = etp * 0.001 * dt
      beta = 1.0
      if (ice .ne. 1) then
        if (esd .lt. etp2) then
          beta = esd / etp2
        endif
      endif
c ----------------------------------------------------------------------
c if etp<0 (downward) then dewfall (=frostfall in this case).
c ----------------------------------------------------------------------
      dew = 0.0
      if (etp .lt. 0.0) then
        dew = -etp * 0.001
      endif
c ----------------------------------------------------------------------
c if precip is falling, calculate heat flux from snow sfc to newly
c accumulating precip.  note that this reflects the flux appropriate for
c the not-yet-updated skin temperature (t1).  assumes temperature of the
c snowfall striking the gound is =sfctmp (lowest model level air temp).
c ----------------------------------------------------------------------
      flx1 = 0.0
      if (snowng) then
        flx1 = cpice * prcp * (t1 - sfctmp)
      else
        if (prcp .gt. 0.0) flx1 = cph2o * prcp * (t1 - sfctmp)
      endif
c ----------------------------------------------------------------------
c calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
c between the snow pack and the soil and on net radiation.
c include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
c fluxes.  flx1 from above, flx2 brought in via commom block rite.
c flx2 reflects freezing rain latent heat flux using t1 calculated in
c penman.
c ----------------------------------------------------------------------
      dsoil = -(0.5 * zsoil(1))
      dtot = snowh + dsoil
      denom = 1.0 + df1 / (dtot * rr * rch)
      t12a = ( (fdown-flx1-flx2-sigma*t24)/rch
     &       + th2 - sfctmp - beta*epsca ) / rr
      t12b = df1 * stc(1) / (dtot * rr * rch)
      t12 = (sfctmp + t12a + t12b) / denom      
c ----------------------------------------------------------------------
c if the 'effective snow-grnd sfc temp' is at or below freezing, no snow
c melt will occur.  set the skin temp to this effective temp.  reduce
c (by sublimation ) or increase (by frost) the depth of the snowpack,
c depending on sign of etp.
c update soil heat flux (ssoil) using new skin temperature (t1)
c since no snowmelt, set accumulated snowmelt to zero, set 'effective'
c precip from snowmelt to zero, set phase-change heat flux from snowmelt
c to zero.
c ----------------------------------------------------------------------
      if (t12 .le. tfreez) then
        t1 = t12
        ssoil = df1 * (t1 - stc(1)) / dtot
        esd = max(0.0, esd-etp2)
        flx3 = 0.0
        ex = 0.0
        snomlt = 0.0
      else
c ----------------------------------------------------------------------
c if the 'effective snow-grnd sfc temp' is above freezing, snow melt
c will occur.  call the snow melt rate,ex and amt, snomlt.  revise the
c effective snow depth.  revise the skin temp because it would have chgd
c due to the latent heat released by the melting. calc the latent heat
c released, flx3. set the effective precip, prcp1 to the snow melt rate,
c ex for use in smflx.  adjustment to t1 to account for snow patches.
c calculate qsat valid at freezing point.  note that esat (saturation
c vapor pressure) value of 6.11e+2 used here is that valid at frzzing
c point.  note that etp from call penman in sflx is ignored here in
c favor of bulk etp over 'open water' at freezing temp.
c update soil heat flux (s) using new skin temperature (t1)
c ----------------------------------------------------------------------
        t1 = tfreez * sncovr + t12 * (1.0 - sncovr)
        qsat = (0.622*6.11e2)/(sfcprs-0.378*6.11e2)
        etp = rch*(qsat-q2)/cp
        etp2 = etp*0.001*dt
        beta = 1.0
        ssoil = df1 * (t1 - stc(1)) / dtot
c ----------------------------------------------------------------------
c if potential evap (sublimation) greater than depth of snowpack.
c beta<1
c snowpack has sublimated away, set depth to zero.
c ----------------------------------------------------------------------
        if (esd .le. etp2) then
          beta = esd / etp2
          esd = 0.0
          ex = 0.0
          snomlt = 0.0
        else
c ----------------------------------------------------------------------
c potential evap (sublimation) less than depth of snowpack, retain
c   beta=1.
c snowpack (esd) reduced by potential evap rate
c etp3 (convert to flux)
c ----------------------------------------------------------------------
          esd = esd-etp2
          etp3 = etp*lsubc
          seh = rch*(t1-th2)
          t14 = t1*t1
          t14 = t14*t14
          flx3 = fdown - flx1 - flx2 - sigma*t14 - ssoil - seh - etp3
          if (flx3 .le .0.0) flx3 = 0.0
          ex = flx3*0.001/lsubf
c ----------------------------------------------------------------------
c snowmelt reduction depending on snow cover
c if snow cover less than 5% no snowmelt reduction
c ***note:  does 'if' below fail to match the melt water with the melt
c           energy?
c ----------------------------------------------------------------------
          if (sncovr .gt. 0.05) ex = ex * sncovr
          snomlt = ex * dt
c ----------------------------------------------------------------------
c esdmin represents a snowpack depth threshold value below which we
c choose not to retain any snowpack, and instead include it in snowmelt.
c ----------------------------------------------------------------------
          if (esd-snomlt .ge. esdmin) then
            esd = esd - snomlt
          else
c ----------------------------------------------------------------------
c snowmelt exceeds snow depth
c ----------------------------------------------------------------------
            ex = esd/dt
            flx3 = ex*1000.0*lsubf
            snomlt = esd
            esd = 0.0
          endif
c ----------------------------------------------------------------------
c end of 'esd .le. etp2' if-block
c ----------------------------------------------------------------------
        endif
        prcp1 = prcp1 + ex
c ----------------------------------------------------------------------
c end of 't12 .le. tfreez' if-block
c ----------------------------------------------------------------------
      endif
c ----------------------------------------------------------------------
c final beta now in hand, so compute evaporation.  evap equals etp
c unless beta<1.
c ----------------------------------------------------------------------
      eta = beta*etp
c ----------------------------------------------------------------------
c set the effective potnl evapotransp (etp1) to zero since this is snow
c case, so surface evap not calculated from edir, ec, or ett in smflx
c (below).
c if seaice (ice=1) skip call to smflx.
c smflx returns updated soil moisture values.  in this, the snow pack
c case, eta1 is not used in calculation of evap.
c ----------------------------------------------------------------------
      etp1 = 0.0
      if (ice .ne. 1) then
        call smflx (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,
     &              sh2o,slope,kdt,frzfact,
     &              smcmax,bexp,pc,smcwlt,dksat,dwsat,
     &              smcref,shdfac,cmcmax,
     &              smcdry,cfactr,runoff1,runoff2,runoff3,
     &              edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,
     &              drip)
      endif
c ----------------------------------------------------------------------
c before call shflx in this snowpack case, set zz1 and yy arguments to
c special values that ensure that ground heat flux calculated in shflx
c matches that already computer for below the snowpack, thus the sfc
c heat flux to be computed in shflx will effectively be the flux at the
c snow top surface.  t11 is a dummy arguement so we will not use the
c skin temp value as revised by shflx.
c ----------------------------------------------------------------------
      zz1 = 1.0
      yy = stc(1)-0.5*ssoil*zsoil(1)*zz1/df1
      t11 = t1
c ----------------------------------------------------------------------
c shflx will calc/update the soil temps.  note:  the sub-sfc heat flux 
c (ssoil1) and the skin temp (t11) output from this shflx call are not
c used  in any subsequent calculations. rather, they are dummy variables
c here in the snopac case, since the skin temp and sub-sfc heat flux are
c updated instead near the beginning of the call to snopac.
c ----------------------------------------------------------------------
      call shflx (ssoil1,stc,smc,smcmax,nsoil,t11,dt,yy,zz1,zsoil,
     &            tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,
     &            quartz,csoil)
c ----------------------------------------------------------------------
c snow depth and density adjustment based on snow compaction.  yy is
c assumed to be the soil temperture at the top of the soil column.
c ----------------------------------------------------------------------
      if (esd .gt. 0.) then
        call snowpack (esd,dt,snowh,sndens,t1,yy)
      else
        esd = 0.
        snowh = 0.
        sndens = 0.
        sncond = 1.
      endif
c ----------------------------------------------------------------------
c end subroutine snopac
c ----------------------------------------------------------------------
      return
      end
