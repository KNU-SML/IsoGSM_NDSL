#define CAL cal
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#define SBC sbc
#define T0C t0c
#include <define.h>
      subroutine sfc_noah(ims2,imx2,kmx,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
     &               drain,evap,hflx,ep,wind,
     &               snwdph, slc, snoalb, slptype, shdmin, shdmax,
#ifndef NOAHYDRO
     &               snowfl, runoff, precip, srflag)
#else
     &               snowfl, runoff, precip, srflag, ec)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
c
c ca is the von karman constant
c
#include <paramodel.h>
      parameter (charnock=.014,ca=.4)
      parameter (rd=rd_,rv=rv_,cp=cp_,g=g_,sigma=sbc_)
      parameter (eps=rd/rv,hvap=hvap_,hfus=hfus_)
      parameter (rvrdm1=rvrdm1_,t0c=t0c_,epsm1=eps-1.)
      parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
      parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)
      parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
      parameter (bb2=-.1954,cc2=.009999)
      parameter (elocp=hvap/cp,dfsnow=.31,ch2o=4.2e6,csoil=1.26e6)
      parameter (scanop=.5,cfactr=.5,zbot=-3.,tgice=271.2)
      parameter (cice=1880.*917.,topt=298.)
      parameter (rhoh2o=1000.,convrad=cal_*1.e4/60.)
      parameter (ctfil1=.5,ctfil2=1.-ctfil1)
      parameter (rnu=1.51e-5,arnu=.135*rnu)
      integer   soiltyp, vegtype, slptype
      integer i, k
      real kt1, kt2, ktsoil
c
c passing array
c
      dimension ps(imx2),t1(imx2),q1(imx2)
      dimension sheleg(imx2),snowmt(imx2),snowev(imx2)
      dimension cm(imx2),ch(imx2)
      dimension tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
      dimension smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
      dimension z0rl(imx2),plantr(imx2),soiltyp(imx2),gflux(imx2)
      dimension slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
      dimension drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
      dimension evap(imx2),hflx(imx2),ep(imx2)
      dimension wind(imx2)
      dimension dlwflx(imx2)
      dimension vegtype(imx2)
      dimension runoff(imx2)
      dimension precip(imx2), srflag(imx2)
      dimension snowfl(imx2)
#ifdef NOAHYDRO
      dimension sn_new(imx2)
#endif
      dimension snoalb(imx2),slptype(imx2),shdmin(imx2),shdmax(imx2)
      dimension snwdph(imx2), slc(imx2, kmx)
c
c local array
c
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      dimension psurf(ILOTS),tsurf(ILOTS),
     *          q0(ILOTS)
      dimension theta1(ILOTS),tv1(ILOTS),
     1           rho(ILOTS)
      dimension qs1(ILOTS),qss(ILOTS),
     1          rch(ILOTS),slwd(ILOTS)
c
      integer   ice, couple, nroot
      dimension sldpth(ILOTS, kmx)
      dimension ice(ILOTS),    zlvl(ILOTS),
     *          radflx(ILOTS)
     &,         sfcprs(ILOTS), prcp(ILOTS), 
     *          dqsdt2(ILOTS)
     &,         snowh(ILOTS), sneqv(ILOTS), 
     *          albedo(ILOTS)
     &,         cmx(ILOTS),     chx(ILOTS),
     *          cmc(ILOTS)
     &,         runoff1(ILOTS),runoff2(ILOTS),
     1          runoff3(ILOTS)
      dimension sldpth1d(kmx), smc1d(kmx)
     &,         slc1d(kmx),    stc1d(kmx)
     &,         et(kmx)
      dimension ptu(ILOTS), alb(ILOTS)
     &,         ec(ILOTS), edir(ILOTS), 
     *          ett(ILOTS)
     &,         esnow(ILOTS), drip(ILOTS),
     *          dew(ILOTS)
     &,         flx1(ILOTS), flx2(ILOTS),
     *          flx3(ILOTS)
     &,         sncovr(ILOTS), rc(ILOTS),
     *          pc(ILOTS)
     &,         rsmin(ILOTS), xlai(ILOTS),
     *          rcs(ILOTS)
     &,         rct(ILOTS), rcq(ILOTS)
     &,         rcsoil(ILOTS),soilw(ILOTS),
     *          soilm(ILOTS)
     &,         smcwlt(ILOTS),smcdry(ILOTS),
     1          smcref(ILOTS)
     &,         smcmax(ILOTS),nroot(ILOTS)
c    
c  declare parameters needed for deriving noah variables
c  source: ETA's surface.f
c
       real a2,a3,a4,a23m4
       parameter(a2=17.2693882,a3=273.16,a4=35.86,a23m4=a2*(a3-a4)) 

c
c  surface energy/water balance over land and seaice
c
c  initialization
c
      latd = 23
      lond = 100
      delt2 = delt * 2.
      im = ims2
      km = kmx
      do i=1, im                                           
       if(slimsk(i).eq.1.) then       
         if(srflag(i) .eq. 1.) then
           snowfl(i)=precip(i)
         else
           snowfl(i)=0.
         endif
       endif
      enddo
c
c  initialize variables. all units are supposedly m.k.s. unless specifie
c  psurf is in pascals
c  wind is wind speed, theta1 is adiabatic surface temp from level 1
c  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at
c  surface
c  convert slrad to the civilized unit from langley minute-1 k-4
c  surface roughness length is converted to m from cm
c
       do i = 1, im
       if(slimsk(i).eq.1.) then       
        psurf(i) = 1000. * ps(i)
        q0(i) = max(q1(i),1.e-8)
        tsurf(i) = tskin(i)
        theta1(i) = t1(i) / slk1
        tv1(i) = t1(i) * (1. + rvrdm1 * q0(i))
        rho(i) = (sl1 * psurf(i)) / (rd * tv1(i))
#ifdef ICE
        qs1(i) = 1000. * fpvs(t1(i))
#else
        qs1(i) = 1000. * fpvs0(t1(i))
#endif
        qs1(i) = eps * qs1(i) / (sl1 * psurf(i) + epsm1 * qs1(i))
        qs1(i) = max(qs1(i), 1.e-8)
        q0(i) = min(qs1(i),q0(i))
#ifdef ICE
        qss(i) = 1000. * fpvs(tsurf(i))
#else
        qss(i) = 1000. * fpvs0(tsurf(i))
#endif
        qss(i) = eps * qss(i) / (psurf(i) + epsm1 * qss(i))
       endif
       enddo
c
c  rcp = rho cp ch v
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
        rch(i) = rho(i) * cp * ch(i) * wind(i)
       endif
      enddo
c
c      gsm/rsm -> noah: prepare variables to run noah lsm
c
c 1. configuration information (c):
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
c
c  1.1   ice        sea-ice flag  (=1: sea-ice, =0: land)
c
       if(slimsk(i).eq.2.) then       
         ice(i) =1 
        else 
         ice(i) =0 
       endif
c
c  1.2   dt         timestep (sec) (dt should not exceed 3600 secs)
c        dt = delt
c  1.3   zlvl       height (m) above ground of atmospheric forcing variables
c
       zlvl(i) = -rd * tv1(i) * log(sl1) / g
c
c  1.4   nsoil      number of soil layers (at least 2)
c
c        nsoil = km 
c  1.5   sldpth     the thickness of each soil layer (m)
c
c  noah:zsoil(1) = -sldpth(1)
c
       sldpth(i,1) = -zsoil(i,1)
       do k = 2,km
c
c  noah: zsoil(k) = -sldpth(k)+zsoil(k-1)
c
       sldpth(i,k) = zsoil(i,k-1) - zsoil(i,k)
       enddo
       endif
      enddo
c 
c 2. forcing data (f):
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
c
c  2.1   lwdn       lw downward radiation (w m-2; positive, not net longwave)
c
c  2.2   soldn      solar downward radiation (w m-2; positive, not net solar)
c
c      lwdn = dlwflx
c
       slwd(i) = slrad(i) * convrad
       radflx(i) = -1. * slwd(i)
c
c  2.3   sfcprs     pressure at height zlvl above ground (pascals)
c
       sfcprs(i) = psurf(i) * sl1
c
c  2.4   prcp       precip rate (kg m-2 s-1)
c
       prcp(i) = rhoh2o * precip(i) / delt
c
c  2.5   sfctmp     air temperature (k) at height zlvl above ground
c        sfctmp =  t1 (temp at lowest model layer, in k)
c  2.6   th2        air potential temperature (k) at height zlvl above ground
c        th2 =  theta1 (see above for theta1 calculations)
c  2.7   q2         mixing ratio at height zlvl above ground (kg kg-1)
c        q2 =  q0 
c
       endif
      enddo

c 
c 3. other forcing (input) data (i):
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
c
c  3.1   sfcspd     wind speed (m s-1) at height zlvl above ground
c        sfcspd = wind (wind speed at lowest model layer, in m/sec)
c  3.2   q2sat      sat mixing ratio at height zlvl above ground (kg kg-1)
c        q2sat =  qs1
c  3.3   dqsdt2     slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1)
c      -eta:  dqsdt(i,j)=qlms(i,j)*a23m4/(tlm(i,j)-a4)**2
c      -eta vs mrfx: qlms(i,j) = q2sat; tlm(i,j)=t1
c
        dqsdt2(i) = qs1(i) * a23m4/(t1(i)-a4)**2
       endif
      enddo
c
c 4. canopy/soil characteristics (s):
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
c
c  4.1   vegtyp     vegetation type (integer index)
c      
c  4.2   soiltyp    soil type (integer index)
c     
c  4.3   slopetyp   class of sfc slope (integer index)
c    
       if(slptype(i) .gt.  9) slptype(i) = 9
c
c  4.4   shdfac     areal fractional coverage of green vegetation (0.0-1.0)
c        shdfac = sigmaf
c  4.5   shdmin     minimum areal fractional coverage of green vegetation
c                (fraction= 0.0-1.0) <= shdfac
c  4.6   ptu        photo thermal unit (plant phenology for annuals/crops)
c                 (not yet used)
c
       ptu(i) = 0.   ! arbitrary value (for future use)
c
c  4.7   alb        backround snow-free surface albedo (fraction)
c      mrf uses 4-type albedo; albedo determined in sflx is not used
c      modify sflx ---> comment out albedo determination
c
       alb(i) = 0.2 
c  4.8   snoalb     upper bound on maximum albedo over deep snow
c    
c  4.9   tbot       bottom soil temperature (local yearly-mean sfc air temp)
c        tbot = tg3 
c
       endif
      enddo
c
c 5. history (state) variables (h):
c
      do i = 1, im
       if(slimsk(i).eq.1.) then       
c
c  5.1  cmc         canopy moisture content (m)
c
        cmc(i) = canopy(i)/rhoh2o  !.. convert from mm to m
c
c  5.2  t1          ground/canopy/snowpack) effective skin temperature (k)
c       t1 = tskin 
c  5.3  stc(nsoil)  soil temp (k)
c       stc(nsoil) = stc(km)
c  5.4  smc(nsoil)  total soil moisture content (volumetric fraction)
c       smc(nsoil) = smc(km)
c  5.5  sh2o(nsoil) unfrozen soil moisture content (volumetric fraction)
c       sh2o(nsoil) = slc(km)
c  5.6  snowh       actual snow depth (m)
c
      snowh(i) = snwdph(i) / rhoh2o     !.. convert from mm to m
c
c  5.7  sneqv       liquid water-equivalent snow depth (m)
c
      sneqv(i) = sheleg(i) / rhoh2o     !.. convert from mm to m
c
c  5.8  albedo      surface albedo including snow effect (unitless fraction)
c     here, albedo is set to 0.2
c     once passed into sflx, it is used to estimate soldn
c
      albedo(i) = 0.2
c
c  5.9  ch          surface exchange coefficient for heat and moisture
c                (m s-1); note: ch is technically a conductance
c     ch and cm are computed in parent model -> input variables
c
      chx(i) = ch(i) * wind(i)            ! compute conductance
c
c  5.10 cm          surface exchange coefficient for momentum (m s-1); note:
c                cm is technically a conductance
c
      cmx(i) = cm(i) * wind(i)            ! compute conductance
      endif
      enddo
c
c    call noah lsm
c
c two new arguments are added in sflx routine
c couple: flag indicating whether to run noah in couple or uncouple mode
c srflag: snow-rain flag (1=snow, 0=rain)
c
        couple = 1            !!<---  run noah lsm in 'couple' mode
        do i = 1, im
        if(slimsk(i) .eq. 1.) then
c
c  from global to local array
c
        do k=1, km
         smc1d(k) = smc(i,k)
         stc1d(k) = stc(i,k)
         slc1d(k) = slc(i,k)
         sldpth1d(k) = sldpth(i,k)
        enddo
        call  sflx (
     &  couple,srflag(i),
     &  ice(i),delt,zlvl(i),km,sldpth1d,
     &  dlwflx(i),radflx(i),sfcprs(i),prcp(i),t1(i),
     &  q0(i),wind(i),
     &  theta1(i),qs1(i),dqsdt2(i),
     &  vegtype(i),soiltyp(i),slptype(i),sigmaf(i),
     &  shdmin(i),ptu(i),alb(i),snoalb(i),tg3(i),
     &  cmc(i),tskin(i),stc1d,smc1d,slc1d,
     &  snowh(i),sneqv(i),albedo(i),chx(i),cmx(i),
     &  evap(i),hflx(i),
     &  ec(i),edir(i),et,ett(i),
     &  esnow(i),drip(i),dew(i),
     &  dm(i), ep(i), gflux(i),
     &  flx1(i),flx2(i),flx3(i),
     &  snowmt(i),sncovr(i),
     &  runoff1(i),runoff2(i),runoff3(i),
     &  rc(i),pc(i),rsmin(i),xlai(i),
     &  rcs(i),rct(i),rcq(i),rcsoil(i),
     &  soilw(i),soilm(i),
#ifndef NOAHYDRO
     &  smcwlt(i),smcdry(i),smcref(i),smcmax(i),nroot(i))
#else
     &  sn_new(i),
     &  smcwlt(i),smcdry(i),smcref(i),smcmax(i),nroot(i))
c
c acr - snowfall must be the same in sflx as in sfc_noah
c
        snowfl(i) = sn_new(i)
#endif
c
c  from local to global array
c
        do k=1, km
         smc(i,k) = smc1d(k)
         stc(i,k) = stc1d(k)
         slc(i,k) = slc1d(k)
        enddo
#ifdef DBG
        print *, 'land diag fields',
     &  lat, i, evap(i),hflx(i),
     &  ec(i),edir(i),et,ett(i),
     &  esnow(i),drip(i),dew(i),
     &  dm(i), ep(i), gflux(i),
     &  flx1(i),flx2(i),flx3(i),
     &  snowmt(i),sncovr(i),
     &  runoff1(i),runoff2(i),runoff3(i)
#endif
        endif 
        enddo
c
c   noah -> gsm/rsm: prepare variables for return to parent model
c
c 6. output (o):
c
c   return the following output fields to parent model
c   eta        actual latent heat flux (w m-2: positive, if upward from sfc)
c   eta = evap
c   sheat      sensible heat flux (w m-2: positive, if upward from sfc)
c   sheat = hflx
c   beta       ratio of actual/potential evap (dimensionless)
c   beta = dm
c   etp        potential evaporation (w m-2)
c   etp = ep
c   ssoil      soil heat flux (w m-2: negative if downward from surface)
c   ssoil = gflux
c   runoff1    surface runoff (m s-1), not infiltrating the surface
c   runoff1 = runoff (in mm s-1) / 1000.
c   runoff2    subsurface runoff (m s-1), drainage out bottom
c   runoff2 = drain (in mm s-1) / 1000.
c
      do i = 1, im
      if(slimsk(i) .eq. 1.) then
c
c  runoff1, runfoff2 (m s-1) -> runoff, drain (mm s-1)
c
        runoff(i) = runoff1(i) * rhoh2o
        drain(i)  = runoff2(i) * rhoh2o
c
c  cmc, snowh, sneqv (m) -> canopy, snwdph, sheleg (mm)
c
        canopy(i) = cmc(i) * rhoh2o
        snwdph(i) = snowh(i) * rhoh2o
        sheleg(i) = sneqv(i) * rhoh2o
#ifdef NOAHYDRO
c
c  transfer snow sublimation (kg/m^2/s)
c
        snowev(i) = esnow(i)
#endif
      endif
      enddo
c
c  compute qsurf and dm
c
      do i = 1, im
      if(slimsk(i) .eq. 1.) then
        qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
        dm(i) = 1.
      endif
      enddo
c
      return
      end
