#include <define.h>
        subroutine sflx(
     &  couple,srflag,
     &  ice,dt,zlvl,nsoil,sldpth,
     &  lwdn,radflx,sfcprs,prcp,sfctmp,q2,sfcspd,
     &  th2,q2sat,dqsdt2,
     &  vegtyp,soiltyp,slopetyp,shdfac,shdmin,ptu,alb,snoalb,tbot,
     &  cmc,t1,stc,smc,sh2o,snowh,sneqv,albedo,ch,cm,
     &  eta,sheat,
c ----------------------------------------------------------------------
c outputs, diagnostics, parameters below generally not necessary when
c coupled with e.g. a nwp model (such as the noaa/nws/ncep mesoscale eta
c model).  other applications may require different output variables. 
c ----------------------------------------------------------------------
     &  ec,edir,et,ett,esnow,drip,dew,
     &  beta,etp,ssoil,
     &  flx1,flx2,flx3,
     &  snomlt,sncovr,
     &  runoff1,runoff2,runoff3,
     &  rc,pc,rsmin,xlai,rcs,rct,rcq,rcsoil,
     &  soilw,soilm,
#ifdef NOAHYDRO
     &  sn_new,
#endif
     &  smcwlt,smcdry,smcref,smcmax,nroot)
c
c cheng-hsuan lu (2002)
c (1) modify calling argument
c     1. add 'couple' (=1: coupled, =0: decoupled)
c     2. add 'srflag' for snow-rain detection
c     3. replace soldn by radflx
c (2) modify source code 
c     1. if couple=1, radflx=fdown 
c        if couple=0, radflx=soldn
c     2. move the code fragment where albedo is determined (alcalc)
c        to the section where sfcdif is called
c     3. invoke the section of alcalc and sfcdif only if couple=0
c     4. use srflag (instead of sfctmp) to determine snow/rain
c     5. change czil from 0.2 to 0.1
c     6. pass smc (instead of slc) to devap
c
      implicit none
c ----------------------------------------------------------------------
c subroutine sflx - version 2.5 - 18 october 2001
c ----------------------------------------------------------------------
c sub-driver for "noah/osu lsm" family of physics subroutines for a
c soil/veg/snowpack land-surface model to update soil moisture, soil
c ice, soil temperature, skin temperature, snowpack water content,
c snowdepth, and all terms of the surface energy balance and surface
c water balance (excluding input atmospheric forcings of downward
c radiation and precip)
c ----------------------------------------------------------------------
c sflx argument list key:
c ----------------------------------------------------------------------
c  c  configuration information
c  f  forcing data
c  i  other (input) forcing data
c  s  surface characteristics
c  h  history (state) variables
c  o  output variables
c  d  diagnostic output
c ----------------------------------------------------------------------
c 1. configuration information (c):
c ----------------------------------------------------------------------
c   couple     couple-uncouple flag  (=1: coupled, =0: decoupled)  **clu_rev**
c   srflag     flag for snow-rain detection)        **clu_rev**
c   ice	       sea-ice flag  (=1: sea-ice, =0: land)
c   dt	       timestep (sec) (dt should not exceed 3600 secs, recommend
c                1800 secs or less)
c   zlvl       height (m) above ground of atmospheric forcing variables
c   nsoil      number of soil layers (at least 2, and not greater than
c                parameter nsold set below)
c   sldpth     the thickness of each soil layer (m)
c ----------------------------------------------------------------------
c 2. forcing data (f):
c ----------------------------------------------------------------------
c   lwdn       lw downward radiation (w m-2; positive, not net longwave)
c   radflx     radiation flux (soldn or fdown)                    **clu_rev**
c   soldn      solar downward radiation (w m-2; positive, not net solar)
c   sfcprs     pressure at height zlvl above ground (pascals)
c   prcp       precip rate (kg m-2 s-1) (note, this is a rate)
c   sfctmp     air temperature (k) at height zlvl above ground
c   th2        air potential temperature (k) at height zlvl above ground
c   q2         mixing ratio at height zlvl above ground (kg kg-1)
c ----------------------------------------------------------------------
c 3. other forcing (input) data (i):
c ----------------------------------------------------------------------
c   sfcspd     wind speed (m s-1) at height zlvl above ground
c   q2sat      sat mixing ratio at height zlvl above ground (kg kg-1)
c   dqsdt2     slope of sat specific humidity curve at t=sfctmp
c                (kg kg-1 k-1)
c ----------------------------------------------------------------------
c 4. canopy/soil characteristics (s):
c ----------------------------------------------------------------------
c   vegtyp     vegetation type (integer index)
c   soiltyp    soil type (integer index)
c   slopetyp   class of sfc slope (integer index)
c   shdfac     areal fractional coverage of green vegetation
c                (fraction= 0.0-1.0)
c   shdmin     minimum areal fractional coverage of green vegetation
c                (fraction= 0.0-1.0) <= shdfac
c   ptu        photo thermal unit (plant phenology for annuals/crops)
c                (not yet used, but passed to redprm for future use in
c                veg parms)
c   alb        backround snow-free surface albedo (fraction), for julian
c                day of year (usually from temporal interpolation of
c                monthly mean values' calling prog may or may not
c                include diurnal sun angle effect)
c   snoalb     upper bound on maximum albedo over deep snow (e.g. from
c                robinson and kukla, 1985, j. clim. & appl. meteor.)
c   tbot       bottom soil temperature (local yearly-mean sfc air
c                temperature)
c ----------------------------------------------------------------------
c 5. history (state) variables (h):
c ----------------------------------------------------------------------
c  cmc         canopy moisture content (m)
c  t1          ground/canopy/snowpack) effective skin temperature (k)
c  stc(nsoil)  soil temp (k)
c  smc(nsoil)  total soil moisture content (volumetric fraction)
c  sh2o(nsoil) unfrozen soil moisture content (volumetric fraction)
c                note: frozen soil moisture = smc - sh2o
c  snowh       actual snow depth (m)
c  sneqv       liquid water-equivalent snow depth (m)
c                note: snow density = sneqv/snowh
c  albedo      surface albedo including snow effect (unitless fraction)
c                =snow-free albedo (alb) when sneqv=0, or
c                =fct(msnoalb,alb,vegtyp,shdfac,shdmin) when sneqv>0
c  ch          surface exchange coefficient for heat and moisture
c                (m s-1); note: ch is technically a conductance since
c                it has been multiplied by wind speed.
c  cm          surface exchange coefficient for momentum (m s-1); note:
c                cm is technically a conductance since it has been
c                multiplied by wind speed.  cm is not needed in sflx
c ----------------------------------------------------------------------
c 6. output (o):
c ----------------------------------------------------------------------
c output variables necessary for a coupled numerical weather prediction
c model, e.g. noaa/nws/ncep mesoscale eta model.  for this application,
c the remaining output/diagnostic/parameter blocks below are not
c necessary.  other applications may require different output variables.
c   eta        actual latent heat flux (kg m-2 s-1: negative, if up from
c	         surface)
c   sheat      sensible heat flux (w m-2: negative, if upward from
c	         surface)
c ----------------------------------------------------------------------
c   ec         canopy water evaporation (m s-1)
c   edir       direct soil evaporation (m s-1)
c   et(nsoil)  plant transpiration from a particular root (soil) layer
c                 (m s-1)
c   ett        total plant transpiration (m s-1)
c   esnow      sublimation from (or deposition to if <0) snowpack
c                (m s-1)
c   drip       through-fall of precip and/or dew in excess of canopy
c                water-holding capacity (m)
c   dew        dewfall (or frostfall for t<273.15) (m)
c ----------------------------------------------------------------------
c   beta       ratio of actual/potential evap (dimensionless)
c   etp        potential evaporation (kg m-2 s-1)
c   ssoil      soil heat flux (w m-2: negative if downward from surface)
c ----------------------------------------------------------------------
c   flx1       precip-snow sfc (w m-2)
c   flx2       freezing rain latent heat flux (w m-2)
c   flx3       phase-change heat flux from snowmelt (w m-2)
c ----------------------------------------------------------------------
c   snomlt     snow melt (m) (water equivalent)
c   sncovr     fractional snow cover (unitless fraction, 0-1)
c ----------------------------------------------------------------------
c   runoff1    surface runoff (m s-1), not infiltrating the surface
c   runoff2    subsurface runoff (m s-1), drainage out bottom of last
c                soil layer
c   runoff3    numerical trunctation in excess of porosity (smcmax)
c                for a given soil layer at the end of a time step
c ----------------------------------------------------------------------
c   rc         canopy resistance (s m-1)
c   pc         plant coefficient (unitless fraction, 0-1) where pc*etp
c                = actual transp
c   xlai       leaf area index (dimensionless)
c   rsmin      minimum canopy resistance (s m-1)
c   rcs        incoming solar rc factor (dimensionless)
c   rct        air temperature rc factor (dimensionless)
c   rcq        atmos vapor pressure deficit rc factor (dimensionless)
c   rcsoil     soil moisture rc factor (dimensionless)
c ----------------------------------------------------------------------
c 7. diagnostic output (d):
c ----------------------------------------------------------------------
c   soilw      available soil moisture in root zone (unitless fraction
c	         between smcwlt and smcmax)
c   soilm      total soil column moisture content (frozen+unfrozen) (m) 
c ----------------------------------------------------------------------
c 8. parameters (p):
c ----------------------------------------------------------------------
c   smcwlt     wilting point (volumetric)
c   smcdry     dry soil moisture threshold where direct evap frm top
c                layer ends (volumetric)
c   smcref     soil moisture threshold where transpiration begins to
c                stress (volumetric)
c   smcmax     porosity, i.e. saturated value of soil moisture
c                (volumetric)
c   nroot      number of root layers, a function of veg type, determined
c              in subroutine redprm.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c ----------------------------------------------------------------------
c declarations - logical
c ----------------------------------------------------------------------
      logical frzgra
      logical saturated
      logical snowng
c ----------------------------------------------------------------------
c declarations - integer
c ----------------------------------------------------------------------
      integer ice
      integer couple                !...**clu_rev**
      integer k
      integer kz
      integer nsoil
      integer nroot
      integer slopetyp
      integer soiltyp
      integer vegtyp
c ----------------------------------------------------------------------
c declarations - real
c ----------------------------------------------------------------------
      real albedo
      real alb
      real bexp
      real beta
      real cfactr
      real ch
      real cm
      real cmc
      real cmcmax
      real cp
      real csnow
      real csoil
      real czil
      real dew
      real df1
      real df1h
      real df1a
      real dksat
      real dt
      real dwsat
      real dqsdt2
      real dsoil
      real dtot
      real drip
      real ec
      real edir
      real esnow
      real et(nsoil)
      real ett
      real frcsno
      real frcsoi
      real epsca
      real eta
      real etp
      real fdown
      real f1
      real flx1
      real flx2
      real flx3
      real fxexp
      real frzx
      real sheat
      real hs
      real kdt
      real lwdn
#ifndef NOAHYDRO
      real lvh2o
#else
      real hvap
#endif
      real pc
      real prcp
      real ptu
      real prcp1
      real psisat
      real q2
      real q2sat
      real quartz
      real r
      real radflx                      ! ..**clu_rev**
      real rch
      real refkdt
      real rr
      real rtdis(nsold)
      real runoff1
      real runoff2
      real rgl
      real runoff3
      real rsmax
      real rc
      real rsmin
      real rcq
      real rcs
      real rcsoil
      real rct
      real rsnow
      real sndens
      real sncond 
      real ssoil
      real sbeta
      real sfcprs
      real sfcspd
      real sfctmp
      real srflag                      ! ..**clu_rev**
      real shdfac
      real shdmin
      real sh2o(nsoil)
      real sldpth(nsoil)
      real smcdry
      real smcmax
      real smcref
      real smcwlt
      real smc(nsoil)
      real sneqv
      real sncovr
      real snowh
      real sn_new
      real slope
      real snup
      real salp
      real snoalb
      real stc(nsoil)
      real snomlt
      real soldn
      real soilm
      real soilw
      real soilwm
      real soilww
      real t1
      real t1v
      real t24
      real t2v
      real tbot
      real th2
      real th2v
      real topt
      real tfreez
      real tsnow
      real xlai
      real zlvl
      real zbot
      real z0
      real zsoil(nsold)
c ----------------------------------------------------------------------
c declarations - parameters
c ----------------------------------------------------------------------
      parameter(tfreez = 273.15)
#ifndef NOAHYDRO
      parameter(lvh2o = 2.501e+6)
#else
c acr - use gsm's latent heat of vaporization
      parameter(hvap=2.500e+6)
#endif
      parameter(r = 287.04)
      parameter(cp = 1004.5)
c ----------------------------------------------------------------------
c   initialization
c ----------------------------------------------------------------------
      runoff1 = 0.0
      runoff2 = 0.0
      runoff3 = 0.0
      snomlt = 0.0
c
	if (ice .eq. 1) then
            print *, 'error:call noah over sea ice land mask '
        endif
c ----------------------------------------------------------------------
c  the variable "ice" is a flag denoting sea-ice case 
c ----------------------------------------------------------------------
      if (ice .eq. 1) then
c ----------------------------------------------------------------------
c sea-ice layers are equal thickness and sum to 3 meters
c ----------------------------------------------------------------------
        do kz = 1,nsoil
          zsoil(kz) = -3.*float(kz)/float(nsoil)
        end do
      else
c ----------------------------------------------------------------------
c calculate depth (negative) below ground from top skin sfc to bottom of
c   each soil layer.  note:  sign of zsoil is negative (denoting below
c   ground)
c ----------------------------------------------------------------------
        zsoil(1) = -sldpth(1)
        do kz = 2,nsoil
          zsoil(kz) = -sldpth(kz)+zsoil(kz-1)
        end do
      endif
c ----------------------------------------------------------------------
c next is crucial call to set the land-surface parameters, including
c soil-type and veg-type dependent parameters.
c ----------------------------------------------------------------------
      call redprm (
     &             vegtyp,soiltyp,slopetyp,
     &      	   cfactr,cmcmax,rsmax,topt,refkdt,kdt,sbeta,
     &      	   shdfac,rsmin,rgl,hs,zbot,frzx,psisat,slope,
     &      	   snup,salp,bexp,dksat,dwsat,smcmax,smcwlt,smcref,
     &      	   smcdry,f1,quartz,fxexp,rtdis,sldpth,zsoil,
     &      	   nroot,nsoil,z0,czil,xlai,csoil,ptu)
c ----------------------------------------------------------------------
c  initialize precipitation logicals.
c ----------------------------------------------------------------------
      snowng = .false.
      frzgra = .false.
c ----------------------------------------------------------------------
c if sea-ice case, assign default water-equiv snow on top
c ----------------------------------------------------------------------
      if (ice .eq. 1) then
        sneqv = 0.01
        snowh = 0.05
      endif
c ----------------------------------------------------------------------
c if input snowpack is nonzero, then compute snow density "sndens" and
c   snow thermal conductivity "sncond" (note that csnow is a function
c   subroutine)
c ----------------------------------------------------------------------
      if (sneqv .eq. 0.0) then
        sndens = 0.0
        snowh = 0.0
        sncond = 1.0
      else
        sndens = sneqv/snowh
c       if(sndens.gt.1.) sndens=1.
        sncond = csnow(sndens) 
      endif
c ----------------------------------------------------------------------
c determine if it's precipitating and what kind of precip it is.
c if it's prcping and the air temp is colder than 0 c, it's snowing!
c if it's prcping and the air temp is warmer than 0 c, but the grnd
c temp is colder than 0 c, freezing rain is presumed to be falling.
c ----------------------------------------------------------------------
      if (prcp .gt. 0.0) then
!clu    if (sfctmp .le. tfreez) then   !** clu_rev **
        if (srflag .eq. 1.) then   
          snowng = .true.
        else
          if (t1 .le. tfreez) frzgra = .true.
        endif
      endif
c ----------------------------------------------------------------------
c if either prcp flag is set, determine new snowfall (converting prcp
c rate from kg m-2 s-1 to a liquid equiv snow depth in meters) and add
c it to the existing snowpack.
c note that since all precip is added to snowpack, no precip infiltrates
c into the soil so that prcp1 is set to zero.
c ----------------------------------------------------------------------
      if ( (snowng) .or. (frzgra) ) then
        sn_new = prcp * dt * 0.001
        sneqv = sneqv + sn_new
        prcp1 = 0.0
c ----------------------------------------------------------------------
c update snow density based on new snowfall, using old and new snow.
c update snow thermal conductivity
c ----------------------------------------------------------------------
        call snow_new (sfctmp,sn_new,snowh,sndens)  
        sncond = csnow (sndens) 
      else
c ----------------------------------------------------------------------
c precip is liquid (rain), hence save in the precip variable that
c later can wholely or partially infiltrate the soil (along with 
c any canopy "drip" added to this later)
#ifdef NOAHYDRO
c acr addition - zero out new snow when it is not snowing
c ----------------------------------------------------------------------
        sn_new = 0
#else
c ----------------------------------------------------------------------
#endif
        prcp1 = prcp
      endif
c
c ** clu_rev **
c only determine snowcover (needed for thermal conductivity calculations)
c the surface albedo modification is done together with sfcdif (optional)
c
c ----------------------------------------------------------------------
c determine snowcover over land.
c ----------------------------------------------------------------------
      if (ice .eq. 0) then
c ----------------------------------------------------------------------
c if snow depth=0, set snow fraction=0
c ----------------------------------------------------------------------
        if (sneqv .eq. 0.0) then
          sncovr = 0.0
!clu      albedo = alb
        else
c ----------------------------------------------------------------------
c determine snow fractional coverage.
c ----------------------------------------------------------------------
          call snfrac (sneqv,snup,salp,snowh,sncovr)
!clu      call alcalc (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
        endif
      else
c ----------------------------------------------------------------------
c set snow cover over sea-ice
c ----------------------------------------------------------------------
        sncovr = 1.0
!clu    albedo = 0.60
      endif
c ----------------------------------------------------------------------
c thermal conductivity for sea-ice case
c ----------------------------------------------------------------------
      if (ice .eq. 1) then
        df1 = 2.2
      else
c ----------------------------------------------------------------------
c next calculate the subsurface heat flux, which first requires
c calculation of the thermal diffusivity.  treatment of the
c latter follows that on pages 148-149 from "heat transfer in 
c cold climates", by v. j. lunardini (published in 1981 
c by van nostrand reinhold co.) i.e. treatment of two contiguous 
c "plane parallel" mediums (namely here the first soil layer 
c and the snowpack layer, if any). this diffusivity treatment 
c behaves well for both zero and nonzero snowpack, including the 
c limit of very thin snowpack.  this treatment also eliminates
c the need to impose an arbitrary upper bound on subsurface 
c heat flux when the snowpack becomes extremely thin.
c ----------------------------------------------------------------------
c first calculate thermal diffusivity of top soil layer, using
c both the frozen and liquid soil moisture, following the 
c soil thermal diffusivity function of peters-lidard et al.
c (1998,jas, vol 55, 1209-1224), which requires the specifying
c the quartz content of the given soil class (see routine redprm)
c ----------------------------------------------------------------------
        call tdfcnd (df1,smc(1),quartz,smcmax,sh2o(1))
c ----------------------------------------------------------------------
c next add subsurface heat flux reduction effect from the 
c overlying green canopy, adapted from section 2.1.2 of 
c peters-lidard et al. (1997, jgr, vol 102(d4))
c ----------------------------------------------------------------------
        df1 = df1 * exp(sbeta*shdfac)
      endif
c ----------------------------------------------------------------------
c finally "plane parallel" snowpack effect following 
c v.j. linardini reference cited above. note that dtot is
c combined depth of snowdepth and thickness of first soil layer
c ----------------------------------------------------------------------
      dsoil = -(0.5 * zsoil(1))
      if (sneqv .eq. 0.) then
        ssoil = df1 * (t1 - stc(1) ) / dsoil
      else
        dtot = snowh + dsoil
        frcsno = snowh/dtot
        frcsoi = dsoil/dtot
c
c 1. harmonic mean (series flow)
c        df1 = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
        df1h = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
c 2. arithmetic mean (parallel flow)
c        df1 = frcsno*sncond + frcsoi*df1
        df1a = frcsno*sncond + frcsoi*df1
c
c 3. geometric mean (intermediate between harmonic and arithmetic mean)
c        df1 = (sncond**frcsno)*(df1**frcsoi)
c test - mbek, 10 jan 2002
c weigh df by snow fraction
c        df1 = df1h*sncovr + df1a*(1.0-sncovr)
c        df1 = df1h*sncovr + df1*(1.0-sncovr)
        df1 = df1a*sncovr + df1*(1.0-sncovr)
c ----------------------------------------------------------------------
c calculate subsurface heat flux, ssoil, from final thermal diffusivity
c of surface mediums, df1 above, and skin temperature and top 
c mid-layer soil temperature
c ----------------------------------------------------------------------
        ssoil = df1 * (t1 - stc(1) ) / dtot
      endif
c ----------------------------------------------------------------------
c determine surface roughness over snowpack using snow condition from
c the previous timestep.
c ----------------------------------------------------------------------
      if (sncovr .gt. 0.) then
        call snowz0 (sncovr,z0)
      endif
c ----------------------------------------------------------------------
c next call routine sfcdif to calculate the sfc exchange coef (ch) for
c heat and moisture.
c
c note !!!
c comment out call sfcdif, if sfcdif already called in calling program
c (such as in coupled atmospheric model).
c
c note !!!
c do not call sfcdif until after above call to redprm, in case
c alternative values of roughness length (z0) and zilintinkevich coef
c (czil) are set there via namelist i/o.
c
c note !!!
c routine sfcdif returns a ch that represents the wind spd times the
c "original" nondimensional "ch" typical in literature.  hence the ch
c returned from sfcdif has units of m/s.  the important companion
c coefficient of ch, carried here as "rch", is the ch from sfcdif times
c air density and parameter "cp".  "rch" is computed in "call penman".
c rch rather than ch is the coeff usually invoked later in eqns.
c
c note !!!
c sfcdif also returns the surface exchange coefficient for momentum, cm,
c also known as the surface drage coefficient, but cm is not used here.
c ----------------------------------------------------------------------
c calc virtual temps and virtual potential temps needed by subroutines
c sfcdif and penman.
c ----------------------------------------------------------------------
      t2v = sfctmp * (1.0 + 0.61 * q2 )
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c ** clu_rev **
c decouple mode
c determine surface albedo 
c compute surface exchange coefficients
c set soldn = radflx
c compute fdown from soldn, albedo, and lwdn
c
c couple mode
c set fdown = radflx
c compute soldn from fdown, albedo, and lwdn
c
      if (couple .eq. 0) then            !......decouple mode
c ----------------------------------------------------------------------
c determine albedo over land.
c ----------------------------------------------------------------------
      if (ice .eq. 0) then
c ----------------------------------------------------------------------
c if snow depth=0, set albedo=snow free albedo.
c ----------------------------------------------------------------------
        if (sneqv .eq. 0.0) then
          albedo = alb
        else
c ----------------------------------------------------------------------
c determine surface albedo modification due to snowdepth state.
c ----------------------------------------------------------------------
         call alcalc (alb,snoalb,shdfac,shdmin,sncovr,tsnow,albedo)
        endif
      else
c ----------------------------------------------------------------------
c set albedo over sea-ice
c ----------------------------------------------------------------------
#ifdef NOAH_JUN05
       albedo = 0.65
#else
       albedo = 0.60
#endif
      endif
c ----------------------------------------------------------------------
c comment out below 2 lines if call sfcdif is commented out, i.e. in the
c coupled model.
c ----------------------------------------------------------------------
c
       t1v = t1 * (1.0 + 0.61 * q2)
       th2v = th2 * (1.0 + 0.61 * q2)
       call sfcdif (zlvl,z0,t1v,th2v,sfcspd,czil,cm,ch)
c ----------------------------------------------------------------------
c calculate total downward radiation (solar plus longwave) needed in
c penman ep subroutine that follows
c ----------------------------------------------------------------------
       soldn = radflx                      !..downward solar radiation  
       fdown = soldn*(1.0-albedo) + lwdn   !..total downward radiation
      else                                 !.....couple mode               
       fdown = radflx                      !..total downward radiation
       soldn = (fdown-lwdn)/(1.0-albedo)   !..dn sw rad 
      endif
c----------------------------------------------------------------------
c call penman subroutine to calculate potential evaporation (etp), and
c other partial products and sums save in common/rite for later
c calculations.
c ----------------------------------------------------------------------
       call penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,t24,ssoil,
     &              q2,q2sat,etp,rch,epsca,rr,snowng,frzgra,
     &              dqsdt2,flx2)
c ----------------------------------------------------------------------
c call canres to calculate the canopy resistance and convert it into pc
c if nonzero greenness fraction
c ----------------------------------------------------------------------
      if (shdfac .gt. 0.) then
c ----------------------------------------------------------------------
c  frozen ground extension: total soil water "smc" was replaced 
c  by unfrozen soil water "sh2o" in call to canres below
c ----------------------------------------------------------------------
        call canres (soldn,ch,sfctmp,q2,sfcprs,sh2o,zsoil,nsoil,
     &               smcwlt,smcref,rsmin,rc,pc,nroot,q2sat,dqsdt2,
     &               topt,rsmax,rgl,hs,xlai,
     &               rcs,rct,rcq,rcsoil)
      endif
c ----------------------------------------------------------------------
c now decide major pathway branch to take depending on whether snowpack
c exists or not:
c ----------------------------------------------------------------------
      esnow = 0.0
      if (sneqv .eq. 0.0) then
        call nopac (etp,eta,prcp,smc,smcmax,smcwlt,
     &     	    smcref,smcdry,cmc,cmcmax,nsoil,dt,shdfac,
     &     	    sbeta,q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,
     &     	    stc,epsca,bexp,pc,rch,rr,cfactr,
     &     	    sh2o,slope,kdt,frzx,psisat,zsoil,
     &     	    dksat,dwsat,tbot,zbot,runoff1,runoff2,
     &     	    runoff3,edir,ec,et,ett,nroot,ice,rtdis,
     &     	    quartz,fxexp,csoil,
     &     	    beta,drip,dew,flx1,flx2,flx3)
      else
        call snopac (etp,eta,prcp,prcp1,snowng,smc,smcmax,smcwlt,
     &               smcref,smcdry,cmc,cmcmax,nsoil,dt,
     &               sbeta,df1,
     &               q2,t1,sfctmp,t24,th2,fdown,f1,ssoil,stc,epsca,
     &               sfcprs,bexp,pc,rch,rr,cfactr,sncovr,sneqv,sndens,
     &               snowh,sh2o,slope,kdt,frzx,psisat,snup,
     &               zsoil,dwsat,dksat,tbot,zbot,shdfac,runoff1,
     &               runoff2,runoff3,edir,ec,et,ett,nroot,snomlt,
     &               ice,rtdis,quartz,fxexp,csoil,
     &               beta,drip,dew,flx1,flx2,flx3)
        esnow = eta
      endif
c ----------------------------------------------------------------------
c   prepare sensible heat (h) for return to parent model
c ----------------------------------------------------------------------
      sheat = -(ch * cp * sfcprs)/(r * t2v) * ( th2 - t1 )
c ----------------------------------------------------------------------
c  convert units and/or sign of total evap (eta), potential evap (etp),
c  subsurface heat flux (s), and runoffs for what parent model expects
c  convert eta from kg m-2 s-1 to w m-2
c ----------------------------------------------------------------------
#ifndef NOAHYDRO
      eta = eta*lvh2o
      etp = etp*lvh2o
#else
      eta = eta*hvap
      etp = etp*hvap
#endif
c ----------------------------------------------------------------------
c convert the sign of soil heat flux so that:
c   ssoil>0: warm the surface  (night time)
c   ssoil<0: cool the surface  (day time)
c ----------------------------------------------------------------------
      ssoil = -1.0*ssoil      
c ----------------------------------------------------------------------
c  convert runoff3 (internal layer runoff from supersat) from m to m s-1
c  and add to subsurface runoff/drainage/baseflow
c ----------------------------------------------------------------------
      runoff3 = runoff3/dt
      runoff2 = runoff2+runoff3
c ----------------------------------------------------------------------
c total column soil moisture in meters (soilm) and root-zone 
c soil moisture availability (fraction) relative to porosity/saturation
c ----------------------------------------------------------------------
      soilm = -1.0*smc(1)*zsoil(1)
      do k = 2,nsoil
        soilm = soilm+smc(k)*(zsoil(k-1)-zsoil(k))
      end do
      soilwm = -1.0*(smcmax-smcwlt)*zsoil(1)
      soilww = -1.0*(smc(1)-smcwlt)*zsoil(1)
      do k = 2,nroot
        soilwm = soilwm+(smcmax-smcwlt)*(zsoil(k-1)-zsoil(k))
        soilww = soilww+(smc(k)-smcwlt)*(zsoil(k-1)-zsoil(k))
      end do
      soilw = soilww/soilwm
c ----------------------------------------------------------------------
c end subroutine sflx
c ----------------------------------------------------------------------
      return
      end
