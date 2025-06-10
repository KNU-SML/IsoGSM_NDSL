      subroutine redprm (
     &                   vegtyp,soiltyp,slopetyp,
     &     		 cfactr,cmcmax,rsmax,topt,refkdt,kdt,sbeta,
     &     		 shdfac,rsmin,rgl,hs,zbot,frzx,psisat,slope,
     &     		 snup,salp,bexp,dksat,dwsat,
     &                   smcmax,smcwlt,smcref,
     &     		 smcdry,f1,quartz,fxexp,rtdis,sldpth,zsoil,
     &     		 nroot,nsoil,z0,czil,lai,csoil,ptu)
#include <paramodel.h>
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine redprm
c ----------------------------------------------------------------------
c internally set (default valuess), or optionally read-in via namelist
c i/o, all soil and vegetation parameters required for the execusion of
c the noah lsm.
c
c optional non-default parameters can be read in, accommodating up to 30
c soil, veg, or slope classes, if the default max number of soil, veg,
c and/or slope types is reset.
c
c future upgrades of routine redprm must expand to incorporate some of
c the empirical parameters of the frozen soil and snowpack physics (such
c as in routines frh2o, snowpack, and snow_new) not yet set in this
c redprm routine, but rather set in lower level subroutines.
c
c set maximum number of soil-, veg-, and slopetyp in data statement.
c ----------------------------------------------------------------------
      integer max_slopetyp
      integer max_soiltyp
      integer max_vegtyp
c
      parameter(max_slopetyp = 30)
      parameter(max_soiltyp = 30)
      parameter(max_vegtyp = 30)
c ----------------------------------------------------------------------
c number of defined soil-, veg-, and slopetyps used.
c ----------------------------------------------------------------------
      integer defined_veg
      integer defined_soil
      integer defined_slope
c
      data defined_veg/13/
      data defined_soil/9/
      data defined_slope/9/
c ----------------------------------------------------------------------
c  set-up soil parameters for given soil type
c  input: soltyp: soil type (integer index)
c  output: soil parameters:
c    maxsmc: max soil moisture content (porosity)
c    refsmc: reference soil moisture (onset of soil moisture
c	     stress in transpiration)
c    wltsmc: wilting pt soil moisture contents
c    drysmc: air dry soil moist content limits
c    satpsi: saturated soil potential
c    satdk:  saturated soil hydraulic conductivity
c    bb:     the 'b' parameter
c    satdw:  saturated soil diffusivity
c    f11:    used to compute soil diffusivity/conductivity
c    quartz:  soil quartz content
c ----------------------------------------------------------------------
c soil types   zobler (1986)	  cosby et al (1984) (quartz cont.(1))
c  1	    coarse	      loamy sand	 (0.82)
c  2	    medium	      silty clay loam	 (0.10)
c  3	    fine	      light clay	 (0.25)
c  4	    coarse-medium     sandy loam	 (0.60)
c  5	    coarse-fine       sandy clay	 (0.52)
c  6	    medium-fine       clay loam 	 (0.35)
c  7	    coarse-med-fine   sandy clay loam	 (0.60)
c  8	    organic	      loam		 (0.40)
c  9	    glacial land ice  loamy sand	 (na using 0.82)
c ----------------------------------------------------------------------
      real bb(max_soiltyp)
      real drysmc(max_soiltyp)
      real f11(max_soiltyp)
      real maxsmc(max_soiltyp)
      real refsmc(max_soiltyp)
      real satpsi(max_soiltyp)
      real satdk(max_soiltyp)
      real satdw(max_soiltyp)
      real wltsmc(max_soiltyp)
      real qtz(max_soiltyp)
c
      real bexp
      real dksat
      real dwsat
      real f1
      real ptu
      real quartz
      real refsmc1
      real smcdry
      real smcmax
      real smcref
      real smcwlt
      real wltsmc1
c ----------------------------------------------------------------------
c soil texture-related arrays.
c ----------------------------------------------------------------------
      data maxsmc/0.421, 0.464, 0.468, 0.434, 0.406, 0.465,
     &  	  0.404, 0.439, 0.421, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      data satpsi/0.04, 0.62, 0.47, 0.14, 0.10, 0.26,
     &  	  0.14, 0.36, 0.04, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
      data satdk /1.41e-5, 0.20e-5, 0.10e-5, 0.52e-5, 0.72e-5,
     &  	  0.25e-5, 0.45e-5, 0.34e-5, 1.41e-5, 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00/
      data bb	 /4.26,  8.72, 11.55, 4.74, 10.73,  8.17,
     &  	  6.77,  5.25,  4.26, 0.00,  0.00,  0.00,
     &  	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00,
     &  	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00,
     &  	  0.00,  0.00,  0.00, 0.00,  0.00,  0.00/
      data qtz   /0.82, 0.10, 0.25, 0.60, 0.52, 0.35,
     &  	  0.60, 0.40, 0.82, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     &  	  0.00, 0.00, 0.00, 0.00, 0.00, 0.00/
c ----------------------------------------------------------------------
c the following 5 parameters are derived later in redprm.f from the soil
c data, and are just given here for reference and to force static
c storage allocation. -dag lohmann, feb. 2001
c ----------------------------------------------------------------------
#ifdef NOAH_JUN05
      data refsmc/0.248, 0.368, 0.398, 0.281, 0.321, 0.361,
     &  	  0.293, 0.301, 0.248, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
#else
      data refsmc/0.283, 0.387, 0.412, 0.312, 0.338, 0.382,
     &            0.315, 0.329, 0.283, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
#endif
      data wltsmc/0.029, 0.119, 0.139, 0.047, 0.100, 0.103,
     &  	  0.069, 0.066, 0.029, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      data drysmc/0.029, 0.119, 0.139, 0.047, 0.100, 0.103,
     &  	  0.069, 0.066, 0.029, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	  0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
      data satdw /5.71e-6, 2.33e-5, 1.16e-5, 7.95e-6, 1.90e-5,
     &  	  1.14e-5, 1.06e-5, 1.46e-5, 5.71e-6, 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00,
     &  	  0.00   , 0.00   , 0.00   , 0.00   , 0.00/
      data f11  /-0.999, -1.116, -2.137, -0.572, -3.201, -1.302,
     &  	 -1.519, -0.329, -0.999,  0.000,  0.000,  0.000,
     &  	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
     &  	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,
     &  	  0.000,  0.000,  0.000,  0.000,  0.000,  0.000/
c ----------------------------------------------------------------------
c set-up vegetation parameters for a given vegetaion type:
c input: vegtyp = vegetation type (integer index)
c ouput: vegetation parameters
c   shdfac: vegetation greenness fraction
c   rsmin:  mimimum stomatal resistance
c   rgl:    parameter used in solar rad term of
c	    canopy resistance function
c   hs:     parameter used in vapor pressure deficit term of
c	    canopy resistance function
c   snup:   threshold snow depth (in water equivalent m) that
c   	    implies 100% snow cover
c ----------------------------------------------------------------------
c ssib vegetation types (dorman and sellers, 1989; jam)
c  1:  broadleaf-evergreen trees  (tropical forest)
c  2:  broadleaf-deciduous trees
c  3:  broadleaf and needleleaf trees (mixed forest)
c  4:  needleleaf-evergreen trees
c  5:  needleleaf-deciduous trees (larch)
c  6:  broadleaf trees with groundcover (savanna)
c  7:  groundcover only (perennial)
c  8:  broadleaf shrubs with perennial groundcover
c  9:  broadleaf shrubs with bare soil
c 10:  dwarf trees and shrubs with groundcover (tundra)
c 11:  bare soil
c 12:  cultivations (the same parameters as for type 7)
c 13:  glacial (the same parameters as for type 11)
c ----------------------------------------------------------------------
      integer nroot
      integer nroot_data(max_vegtyp)
c
      real frzfact
      real hs
      real hstbl(max_vegtyp)
      real lai
      real lai_data(max_vegtyp)
      real psisat
      real rsmin
      real rgl
      real rgltbl(max_vegtyp)
      real rsmtbl(max_vegtyp)
      real shdfac
      real snup
      real snupx(max_vegtyp)
      real z0
      real z0_data(max_vegtyp)
c ----------------------------------------------------------------------
c vegetation class-related arrays
c ----------------------------------------------------------------------
      data nroot_data /4,4,4,4,4,4,3,3,3,2,3,3,2,0,0,
     &                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
#ifdef NOAH_JUN05
      data rsmtbl /300.0, 175.0, 175.0, 300.0, 300.0, 70.0,
     &  	    45.0, 225.0, 225.0, 225.0, 400.0, 45.0,
     &  	   150.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     &  	     0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     &  	     0.0,   0.0,   0.0,   0.0,   0.0,  0.0/
#else
      data rsmtbl /150.0, 100.0, 125.0, 150.0, 100.0, 70.0,
     &              40.0, 300.0, 400.0, 150.0, 400.0, 40.0,
     &             150.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     &               0.0,   0.0,   0.0,   0.0,   0.0,  0.0/
#endif
      data rgltbl /30.0,  30.0,  30.0,  30.0,  30.0,  65.0,
     &  	  100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
     &  	  100.0,   0.0,   0.0,   0.0,	0.0,   0.0,
     &  	    0.0,   0.0,   0.0,   0.0,	0.0,   0.0,
     &  	    0.0,   0.0,   0.0,   0.0,	0.0,   0.0/
      data hstbl /41.69, 54.53, 51.93, 47.35,  47.35, 54.53,
     &  	  36.35, 42.00, 42.00, 42.00,  42.00, 36.35,
     &  	  42.00,  0.00,  0.00,  0.00,	0.00,  0.00,
     &  	   0.00,  0.00,  0.00,  0.00,	0.00,  0.00,
     &  	   0.00,  0.00,  0.00,  0.00,	0.00,  0.00/
#ifdef NOAH_JUN05
      data snupx  /0.040, 0.040, 0.040, 0.040, 0.040, 0.040,
     &  	   0.020, 0.020, 0.020, 0.020, 0.013, 0.020,
     &  	   0.013, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	   0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	   0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
#else
      data snupx  /0.080, 0.080, 0.080, 0.080, 0.080, 0.080,
     &             0.040, 0.040, 0.040, 0.040, 0.025, 0.040,
     &             0.025, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &             0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
#endif
      data z0_data /2.653, 0.826, 0.563, 1.089, 0.854, 0.856,
     &  	    0.035, 0.238, 0.065, 0.076, 0.011, 0.035,
     &  	    0.011, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	    0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &  	    0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
#ifdef NOAH_JUN05
      data lai_data /3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     &  	     3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     &  	     3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &  	     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &  	     0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
#else
      data lai_data /4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
     &               4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
     &               4.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &               0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
#endif
c ----------------------------------------------------------------------
c class parameter 'slopetyp' was included to estimate linear reservoir
c coefficient 'slope' to the baseflow runoff out of the bottom layer.
c lowest class (slopetyp=0) means highest slope parameter = 1.
c definition of slopetyp from 'zobler' slope type:
c slope class  percent slope
c 1	       0-8
c 2	       8-30
c 3	       > 30
c 4	       0-30
c 5	       0-8 & > 30
c 6	       8-30 & > 30
c 7	       0-8, 8-30, > 30
c 9	       glacial ice
c blank        ocean/sea
c ----------------------------------------------------------------------
c note:
c class 9 from 'zobler' file should be replaced by 8 and 'blank' 9
c ----------------------------------------------------------------------
      real slope
      real slope_data(max_slopetyp)
c
      data slope_data /0.1,  0.6, 1.0, 0.35, 0.55, 0.8,
     &  	       0.63, 0.0, 0.0, 0.0,  0.0,  0.0,
     &  	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,
     &  	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0,
     &  	       0.0 , 0.0, 0.0, 0.0,  0.0,  0.0/
c ----------------------------------------------------------------------
c set namelist file name
c ----------------------------------------------------------------------
      character*50 namelist_name
c ----------------------------------------------------------------------
c set universal parameters (not dependent on soil, veg, slope type)
c ----------------------------------------------------------------------
      integer i
      integer nsoil
      integer slopetyp
      integer soiltyp
      integer vegtyp
c
      integer bare
      data bare /11/
c
      logical lparam
      data lparam /.true./
c
      logical lfirst
      data lfirst /.true./
c
c ** clu_rev: retain definition status to quarantee 'one-time' execution
      save  lfirst
c
c ** clu_rev: change the default value of czil from 0.2 to 0.1
c ----------------------------------------------------------------------
c parameter used to calculate roughness length of heat.
c ----------------------------------------------------------------------
      real czil
      real czil_data
!clu  data czil_data /0.2/
#ifdef NOAH_JUN05
      data czil_data /0.075/
#else
      data czil_data /0.1/
#endif
c ----------------------------------------------------------------------
c parameter used to caluculate vegetation effect on soil heat flux.
c ----------------------------------------------------------------------
      real sbeta
      real sbeta_data
      data sbeta_data /-2.0/
c ----------------------------------------------------------------------
c bare soil evaporation exponent used in devap.
c ----------------------------------------------------------------------
      real fxexp
      real fxexp_data
      data fxexp_data /2.0/
c ----------------------------------------------------------------------
c soil heat capacity [j m-3 k-1]
c ----------------------------------------------------------------------
      real csoil
      real csoil_data
c     data csoil_data /1.26e+6/
      data csoil_data /2.00e+6/
c ----------------------------------------------------------------------
c specify snow distribution shape parameter salp - shape parameter of
c distribution function of snow cover. from anderson's data (hydro-17)
c best fit is when salp = 2.6
c ----------------------------------------------------------------------
      real salp
      real salp_data
#ifdef NOAH_JUN05
      data salp_data /4.0/
#else
      data salp_data /2.6/
#endif
c ----------------------------------------------------------------------
c kdt is defined by reference refkdt and dksat; refdk=2.e-6 is the sat.
c dk. value for the soil type 2
c ----------------------------------------------------------------------
      real refdk
      real refdk_data
      data refdk_data /2.0e-6/
c
      real refkdt
      real refkdt_data
      data refkdt_data /3.0/
c
      real frzx
      real kdt
c ----------------------------------------------------------------------
c frozen ground parameter, frzk, definition: ice content threshold above
c which frozen soil is impermeable reference value of this parameter for
c the light clay soil (type=3) frzk = 0.15 m.
c ----------------------------------------------------------------------
      real frzk
      real frzk_data
      data frzk_data /0.15/
c
      real rtdis(nsoil)
      real sldpth(nsoil)
      real zsoil(nsoil)
c ----------------------------------------------------------------------
c set two canopy water parameters.
c ----------------------------------------------------------------------
      real cfactr
      real cfactr_data
      data cfactr_data /0.5/
c
      real cmcmax
      real cmcmax_data
      data cmcmax_data /0.5e-3/
c ----------------------------------------------------------------------
c set max. stomatal resistance.
c ----------------------------------------------------------------------
      real rsmax
      real rsmax_data
      data rsmax_data /5000.0/
c ----------------------------------------------------------------------
c set optimum transpiration air temperature.
c ----------------------------------------------------------------------
      real topt
      real topt_data
      data topt_data /298.0/
c ----------------------------------------------------------------------
c specify depth[m] of lower boundary soil temperature.
c ----------------------------------------------------------------------
      real zbot
      real zbot_data
c      data zbot_data /-3.0/
      data zbot_data /-8.0/
c ----------------------------------------------------------------------
c set two soil moisture wilt, soil moisture reference parameters
c ----------------------------------------------------------------------
      real smlow
      real smlow_data
      data smlow_data /0.5/
c
      real smhigh
      real smhigh_data
#ifdef NOAH_JUN05
      data smhigh_data /6.0/
#else
      data smhigh_data /3.0/
#endif
c ----------------------------------------------------------------------
c namelist definition:
c ----------------------------------------------------------------------
      namelist /soil_veg/ slope_data, rsmtbl, rgltbl, hstbl, snupx,
     &  bb, drysmc, f11, maxsmc, refsmc, satpsi, satdk, satdw,
     &  wltsmc, qtz, lparam, zbot_data, salp_data, cfactr_data,
     &  cmcmax_data, sbeta_data, rsmax_data, topt_data,
     &  refdk_data, frzk_data, bare, defined_veg, defined_soil,
     &  defined_slope, fxexp_data, nroot_data, refkdt_data, z0_data,
     &  czil_data, lai_data, csoil_data
c ----------------------------------------------------------------------
c read namelist file to override default parameters only once.
c namelist_name must be 50 characters or less.
c ----------------------------------------------------------------------
      if (lfirst) then
c ** clu_rev: hardwire the namelist**
c        namelist_name = 'soil_veg_namelist_ver_2.5.1'
c
      lfirst = .false.
c
      if (defined_soil .gt. max_soiltyp) then
        write(*,*) 'warning: defined_soil too large in namelist'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
      if (defined_veg .gt. max_vegtyp) then
        write(*,*) 'warning: defined_veg too large in namelist'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
      if (defined_slope .gt. max_slopetyp) then
        write(*,*) 'warning: defined_slope too large in namelist'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
c
      if(vegtyp.eq.0) then
        print *,'vegtyp.eq.0'
      endif
      if(soiltyp.eq.0) then
        print *,'soiltyp.eq.0'
      endif
      if(slopetyp.eq.0) then
        print *,'slopetyp.eq.0'
      endif
c
        smlow = smlow_data
        smhigh = smhigh_data
c
      do i = 1,defined_soil
      satdw(i)  = bb(i)*satdk(i)*(satpsi(i)/maxsmc(i))
      f11(i) = alog10(satpsi(i)) + bb(i)*alog10(maxsmc(i)) + 2.0
      refsmc1 = maxsmc(i)*(5.79e-9/satdk(i))
     &      **(1.0/(2.0*bb(i)+3.0))
c          refsmc(i) = refsmc1 + (maxsmc(i)-refsmc1) / 3.0
      refsmc(i) = refsmc1 + (maxsmc(i)-refsmc1) / smhigh
          wltsmc1 = maxsmc(i) * (200.0/satpsi(i))**(-1.0/bb(i))
c          wltsmc(i) = wltsmc1 - 0.5 * wltsmc1
      wltsmc(i) = wltsmc1 - smlow * wltsmc1
c ----------------------------------------------------------------------
c current version drysmc values that equate to wltsmc.
c future version could let drysmc be independently set via namelist.
c ----------------------------------------------------------------------
      drysmc(i) = wltsmc(i)
      end do

c ----------------------------------------------------------------------
c end lfirst block
c ----------------------------------------------------------------------
      endif
c
      if (soiltyp .gt. defined_soil) then
        write(*,*) 'warning: too many soil types'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
      if (vegtyp .gt. defined_veg) then
        write(*,*) 'warning: too many veg types'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
      if (slopetyp .gt. defined_slope) then
        write(*,*) 'warning: too many slope types'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
c ----------------------------------------------------------------------
c set-up universal parameters (not dependent on soiltyp, vegtyp or
c slopetyp)
c ----------------------------------------------------------------------
      zbot = zbot_data
      salp = salp_data
      cfactr = cfactr_data
      cmcmax = cmcmax_data
      sbeta = sbeta_data
      rsmax = rsmax_data
      topt = topt_data
      refdk = refdk_data
      frzk = frzk_data
      fxexp = fxexp_data
      refkdt = refkdt_data
      czil = czil_data
      csoil = csoil_data
c ----------------------------------------------------------------------
c  set-up soil parameters
c ----------------------------------------------------------------------
      bexp = bb(soiltyp)
      dksat = satdk(soiltyp)
      dwsat = satdw(soiltyp)
      f1 = f11(soiltyp)
c     frzfact = (smcmax / smcref) * (0.412 / 0.468)
      kdt = refkdt * dksat/refdk
      psisat = satpsi(soiltyp)
      quartz = qtz(soiltyp)
      smcdry = drysmc(soiltyp)
      smcmax = maxsmc(soiltyp)
      smcref = refsmc(soiltyp)
      smcwlt = wltsmc(soiltyp)
      frzfact = (smcmax / smcref) * (0.412 / 0.468)
c ----------------------------------------------------------------------
c to adjust frzk parameter to actual soil type: frzk * frzfact
c ----------------------------------------------------------------------
      frzx = frzk * frzfact
c ----------------------------------------------------------------------
c set-up vegetation parameters
c ----------------------------------------------------------------------
      nroot = nroot_data(vegtyp)
c
clu .... force noah to run on 2 layers
c
      if(nroot .gt. nsoil) nroot = nsoil
c
clu ..................................
c
      snup = snupx(vegtyp)
      rsmin = rsmtbl(vegtyp)
      rgl = rgltbl(vegtyp)
      hs = hstbl(vegtyp)
      z0 = z0_data(vegtyp)
      lai = lai_data(vegtyp)
      if (vegtyp .eq. bare) shdfac = 0.0
c
      if (nroot .gt. nsoil) then
        write(*,*) 'warning: too many root layers'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
c ----------------------------------------------------------------------
c calculate root distribution.  present version assumes uniform
c distribution based on soil layer depths.
c ----------------------------------------------------------------------
      do i = 1,nsoil
      rtdis(i) = -sldpth(i)/zsoil(nroot)
      end do
c ----------------------------------------------------------------------
c  set-up slope parameter
c ----------------------------------------------------------------------
      slope = slope_data(slopetyp)
c ----------------------------------------------------------------------
c end subroutine redprm
c ----------------------------------------------------------------------
      return
      end
