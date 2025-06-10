      subroutine nopac(etp,eta,prcp,smc,smcmax,smcwlt,
     &                 smcref,smcdry,cmc,cmcmax,nsoil,dt,shdfac,
     &                 sbeta,q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,
     &                 stc,epsca,bexp,pc,rch,rr,cfactr, 
     &                 sh2o,slope,kdt,frzfact,psisat,zsoil,
     &                 dksat,dwsat,tbot,zbot,runoff1,runoff2,
     &                 runoff3,edir,ec,et,ett,nroot,ice,rtdis,
     &                 quartz,fxexp,csoil,
     &                 beta,drip,dew,flx1,flx2,flx3)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine nopac
c ----------------------------------------------------------------------
c calculate soil moisture and heat flux values and update soil moisture
c content and soil heat content values for the case when no snow pack is
c present.
c ----------------------------------------------------------------------
      integer ice
      integer nroot
      integer nsoil
c
      real bexp
      real beta
      real cfactr
      real cmc
      real cmcmax
      real cp
      real csoil
      real dew
      real df1
      real dksat
      real drip
      real dt
      real dwsat
      real ec
      real edir
      real epsca
      real eta
      real eta1
      real etp
      real etp1
      real et(nsoil)
      real ett
      real fdown
      real f1
      real fxexp
      real flx1
      real flx2
      real flx3
      real frzfact
      real kdt
      real pc
      real prcp
      real prcp1
      real psisat
      real q2
      real quartz
      real rch
      real rr
      real rtdis(nsoil)
      real runoff1
      real runoff2
      real runoff3
      real ssoil
      real sbeta
      real sfctmp
      real shdfac
      real sh2o(nsoil)
      real sigma
      real slope
      real smc(nsoil)
      real smcdry
      real smcmax
      real smcref
      real smcwlt
      real stc(nsoil)
      real t1
      real t24
      real tbot
      real th2
      real yy
      real yynum
      real zbot
      real zsoil(nsoil)
      real zz1
c
      parameter(cp = 1004.5)
      parameter(sigma = 5.67e-8)
c ----------------------------------------------------------------------
c executable code begins here:
c convert etp from kg m-2 s-1 to ms-1 and initialize dew.
c ----------------------------------------------------------------------
      prcp1 = prcp * 0.001
      etp1 = etp * 0.001
      dew = 0.0
      if (etp .gt. 0.0) then
c ----------------------------------------------------------------------
c convert prcp from 'kg m-2 s-1' to 'm s-1'.
c ----------------------------------------------------------------------
           call smflx (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,
     &                 sh2o,slope,kdt,frzfact,
     &                 smcmax,bexp,pc,smcwlt,dksat,dwsat,
     &                 smcref,shdfac,cmcmax,
     &                 smcdry,cfactr,runoff1,runoff2,runoff3, 
     &                 edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,
     &                 drip)
c ----------------------------------------------------------------------
c       convert modeled evapotranspiration fm  m s-1  to  kg m-2 s-1
c ----------------------------------------------------------------------
        eta = eta1 * 1000.0
      else
c ----------------------------------------------------------------------
c if etp < 0, assume dew forms (transform etp1 into dew and reinitialize
c etp1 to zero).
c ----------------------------------------------------------------------
        dew = -etp1
        etp1 = 0.0
c ----------------------------------------------------------------------
c convert prcp from 'kg m-2 s-1' to 'm s-1' and add dew amount.
c ----------------------------------------------------------------------
        prcp1 = prcp1 + dew
c
      call smflx (eta1,smc,nsoil,cmc,etp1,dt,prcp1,zsoil,
     &            sh2o,slope,kdt,frzfact,
     &            smcmax,bexp,pc,smcwlt,dksat,dwsat,
     &            smcref,shdfac,cmcmax,
     &            smcdry,cfactr,runoff1,runoff2,runoff3, 
     &            edir,ec,et,ett,sfctmp,q2,nroot,rtdis,fxexp,
     &            drip)
c ----------------------------------------------------------------------
c convert modeled evapotranspiration from 'm s-1' to 'kg m-2 s-1'.
c ----------------------------------------------------------------------
        eta = eta1 * 1000.0
      endif
c ----------------------------------------------------------------------
c based on etp and e values, determine beta
c ----------------------------------------------------------------------
      if ( etp .le. 0.0 ) then
        beta = 0.0
        if ( etp .lt. 0.0 ) then
          beta = 1.0
          eta = etp
        endif
      else
        beta = eta / etp
      endif
c ----------------------------------------------------------------------
c get soil thermal diffuxivity/conductivity for top soil lyr,
c calc. adjusted top lyr soil temp and adjusted soil flux, then
c call shflx to compute/update soil heat flux and soil temps.
c ----------------------------------------------------------------------
      call tdfcnd (df1,smc(1),quartz,smcmax,sh2o(1))
c ----------------------------------------------------------------------
c vegetation greenness fraction reduction in subsurface heat flux 
c via reduction factor, which is convenient to apply here to thermal 
c diffusivity that is later used in hrt to compute sub sfc heat flux
c (see additional comments on veg effect sub-sfc heat flx in 
c routine sflx)
c ----------------------------------------------------------------------
      df1 = df1 * exp(sbeta*shdfac)
c ----------------------------------------------------------------------
c compute intermediate terms passed to routine hrt (via routine 
c shflx below) for use in computing subsurface heat flux in hrt
c ----------------------------------------------------------------------
      yynum = fdown - sigma * t24
      yy = sfctmp + (yynum/rch+th2-sfctmp-beta*epsca) / rr
      zz1 = df1 / ( -0.5 * zsoil(1) * rch * rr ) + 1.0
      call shflx (ssoil,stc,smc,smcmax,nsoil,t1,dt,yy,zz1,zsoil,
     &            tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,
     &            quartz,csoil)
c ----------------------------------------------------------------------
c set flx1 and flx3 (snopack phase change heat fluxes) to zero since
c they are not used here in snopac.  flx2 (freezing rain heat flux) was
c similarly initialized in the penman routine.
c ----------------------------------------------------------------------
      flx1 = 0.0
      flx3 = 0.0
c ----------------------------------------------------------------------
c end subroutine nopac
c ----------------------------------------------------------------------
      return
      end
