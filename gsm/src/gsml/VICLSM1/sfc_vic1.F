#include <define.h>
      subroutine sfc_vic1(ims2,imx2,msl,lat,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smcld,stcld,dm,sigmaf,vegtype,canopy,
     &               dlwflx,dswflx,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,sl1,slk1,slimsk,
     &               drain,evap,hflx,ep,wind,
     +               nsl,month,binf,ds,dsm,ws,cef,
     +               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,
     +               wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,
     +               apnld,btnld,gmnld,flaild,vrtld,lstsn,
     +               silz,snwz,csno,rsno,tsfld,tpkld,sfwld,pkwld,
     &               snowfl, runoff, precip)
c
#include <paramodel.h>
c
      integer lat
      real rvrdm1, rd, eps, epsm1, cp, g
      real elocp, sigma, hvap
      parameter (rvrdm1=0.6077338, rd=2.8705e+2,eps=0.622)
      parameter (epsm1=-0.378,cp=1.0046e+3,g=9.8)
      parameter (elocp=2488.55,sigma=5.6730e-8,hvap=2.5000e+6)
c    
c
c  notes: vic1 uses one vegetation cover. we will consider multi vegcov
c         in vic2. to simplify coupling vic1/vic2 to gsm/rsm, we use the
c         following local variables to convert the related variables to
c         the uwvic variables. ecpc/ji
c
c  the following definition is for driving vic
c
* ------------------- input variables ----------------------------------

* -- model basic parameter
      
      integer ims2, imx2
      integer msl         ! number of soil layers
      integer nsl         ! number of soil nodes
      integer msub        ! maximum number of subgrid points
      integer md          ! md (1no/2with) precipiation distribution

      parameter (msub=2)  ! for viclsm1 (=2:veg,noveg)
      parameter (md=1)    ! for viclsm1

      real dtime          ! vic time step (second)
      integer month       ! current month

* -- atmosphere data

      real prc(imx2)      ! precipitation rate (m h2o/time step)
      real pgcm(imx2)     ! atm bottom level pressure (pa)
      real wind(imx2)     ! atm bottom level wind speed (m/s)
      real tgcm(imx2)     ! atm bottom level temperature (kelvin)
      real qgcm(imx2)     ! atm bottom level specific humidity (kg/kg)
      real zgcm(imx2)     ! atm bottom level height above surface (m)
      real flwds(imx2)    ! downward longwave rad onto surface (w/m2)
      real sols(imx2)     ! solar rad onto surface (w/m2)

* -- soil parameters

      real binf(imx2)       ! variable infiltration curve parameter (n/a)
      real ds(imx2)         ! fract of dsm nonlinear baseflow begins (fract)
      real dsm(imx2)        ! maximum velocity of baseflow (mm/day)
      real ws(imx2)         ! fract maxi sm nonlinear baseflow occurs (fract)
      real cef(imx2)        ! exponent used in infiltration curve (n/a)
      real silz(imx2)       ! surface roughness of bare soil (m)
      real snwz(imx2)       ! surface roughness of snowpack (m)

      real expt(msl,imx2)   ! para the vari of ksat with sm (n/a)
      real kst(msl,imx2)    ! saturated hydrologic conductivity (mm/day)
      real dph(msl,imx2)    ! thickness of soil layer (m)
      real bub(msl,imx2)    ! bubbling pressure of soil layer (cm)
      real qrt(msl,imx2)    ! quartz content of soil layer (fraction)
      real bkd(msl,imx2)    ! bulk density of soil layer (kg/m3)
      real sld(msl,imx2)    ! soil density of soil layer (kg/m3)
      real wcr(msl,imx2)    ! sm content at the critical point (mm)
      real wpw(msl,imx2)    ! sm content wilting point (mm)
      real smr(msl,imx2)    ! soil moisture residual moisture (mm)
      real smx(msl,imx2)    ! maximum soil moisture (mm)
      real dphnd(nsl,imx2)  ! thickness of soil node (m)
      real smxnd(nsl,imx2)  ! maximum soil moisture at soil node (m3/m3)
      real expnd(nsl,imx2)  ! para the vari of ksat at soil node (n/a)
      real bubnd(nsl,imx2)  ! bubbling pressure at soil node (cm)
      real alpnd(nsl,imx2)  ! para alpha at soil node
      real betnd(nsl,imx2)  ! para beta at soil node
      real gamnd(nsl,imx2)  ! para gamma at soil node

* -- vegetation parameters

      integer nveg(imx2)      ! number of vegetation type in a grid cell
      integer nvt(msub,imx2)  ! vegetation type number
      real wt(msub,imx2)      ! fraction of grid cell covered by veg
      real rt(msl,msub,imx2)  ! fraction of root in the soil layer
      real flai(msub,imx2)    ! leaf area index

* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real cwt(md,msub,imx2)     ! canopy water (m h2o)
      real csn(md,msub,imx2)     ! canopy snow  (m h2o)
      real smc(msl,md,msub,imx2) ! soil moisture content (mm)
      real sic(msl,md,msub,imx2) ! soil ice content (mm)
      real tnd(nsl,md,msub,imx2) ! temperature at soil nodes (k)
      real swq(md,msub,imx2)     ! snow water equivalent (m h2o)
      real rsn(md,msub,imx2)     ! snow density (kg/m^3)
      real tsf(md,msub,imx2)     ! snow surface temperature (k)
      real tpk(md,msub,imx2)     ! snow pack temperature (k)
      real sfw(md,msub,imx2)     ! surface snow water equivalent (m h2o)
      real pkw(md,msub,imx2)     ! snow pack snow water equivalent (m h2o)
      real fmu(imx2)             ! precipitation fraction
      integer lstsn(imx2)        ! number of time step since last snow
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real sh(imx2)       ! sensible heat flux (w/m**2) [+ to ground]
      real lh(imx2)       ! latent heat flux (w/m**2) [+ to ground]
      real gh(imx2)       ! ground heat flux (w/m**2) [+ to ground]
      real ts(imx2)       ! surface radiative temperature (kelvin)
      real albd(imx2)     ! albedo (fraction)
      real ovflw(imx2)    ! overland flow (mm/time step)
      real bsflw(imx2)    ! base flow (mm/time step)
      real snowmt(imx2)   ! snow melt (m/s)
      real snowev(imx2)   ! evaporation over snow surface (w/m2)
      real snowfl(imx2)   ! snow fall (m/time step)
* ----------------------------------------------------------------------
c
c passing array
c
      integer vegtype

      dimension ps(imx2),t1(imx2),q1(imx2),sheleg(imx2)
      dimension tskin(imx2),qsurf(imx2),dm(imx2),smcld(imx2,msl)
      dimension stcld(imx2,nsl),sigmaf(imx2),vegtype(imx2)
      dimension canopy(imx2),dlwflx(imx2),dswflx(imx2)
      dimension z0rl(imx2),tg3(imx2),gflux(imx2)
      dimension zsoil(imx2,msl),cm(imx2),ch(imx2),slimsk(imx2)
      dimension drain(imx2),evap(imx2),hflx(imx2),ep(imx2)
      dimension expld(imx2,msl),dphld(imx2,msl)
      dimension bubld(imx2,msl),qrtld(imx2,msl),bkdld(imx2,msl)
      dimension sldld(imx2,msl),wcrld(imx2,msl),wpwld(imx2,msl)
      dimension smrld(imx2,msl),smxld(imx2,msl),sicld(imx2,msl)
      dimension dpnld(imx2,nsl),sxnld(imx2,nsl),epnld(imx2,nsl)
      dimension bbnld(imx2,nsl),apnld(imx2,nsl),btnld(imx2,nsl)
      dimension gmnld(imx2,nsl),flaild(imx2),vrtld(imx2,msl)
      dimension csno(imx2),rsno(imx2)
      dimension tsfld(imx2),tpkld(imx2),sfwld(imx2),pkwld(imx2)
      dimension sfall(imx2),runoff(imx2),precip(imx2)
      real      kstld(imx2,msl)
c
c local array
c
      dimension psurf(imx2),q0(imx2),theta1(imx2)
      dimension tv1(imx2),rho(imx2),qs1(imx2),rch(imx2)
c
c  surface energy/water balance over land and seaice
c
c  initialization
c

      im = ims2
c
c  initialize variables. all units are supposedly m.k.s. unless specifie
c  psurf is in pascals
c  wind is wind speed, theta1 is adiabatic surface temp from level 1
c  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at surface
c  surface roughness length is converted to m from cm
c
#ifdef DBG
      print*,' enter sfc_vic1 '
#endif
c
      do i = 1, im
         if(slimsk(i).eq.1.) then       
            psurf(i) = 1000. * ps(i)
            q0(i) = max(q1(i),1.e-8)
            theta1(i) = t1(i) / slk1
            tv1(i) = t1(i) * (1. + rvrdm1 * q0(i))
            rho(i) = (sl1 * psurf(i)) / (rd * tv1(i))
#ifdef ICE
            qs1(i) = 1000. * fpvs(t1(i))
#else 
            qs1(i) = 1000. * fpvs0(t1(i))
#endif
            qs1(i) = eps*qs1(i)/(sl1*psurf(i)+epsm1*qs1(i))
            qs1(i) = max(qs1(i), 1.e-8)
            q0(i) = min(qs1(i),q0(i))
c
c  rcp = rho cp ch v
c
            rch(i) = rho(i) * cp * ch(i) * wind(i)
         endif
      enddo
c
c     gsm/rsm -> vic: prepare variables to run vic lsm
c
c 1. model basic parameters
c
c 1.1  msl      number of soil layers
c 1.2  nsl      number of soil nodes
c 1.3  msub     number of subgrids
c 1.4  md       number of wet and dry 
c
c 1.5  dtime    vic lsm time step (sec)
c
      dtime = delt
c
c 1.6  mm       current month
c      
      mm = month
c
c
c 2. atmosphere forcings
c
      do i = 1, im
         if(slimsk(i).eq.1.) then       
c
c 2.1  prc      precipitation rate (m h2o/time step)
c     
            prc(i) = precip(i)
c
c 2.2  pgcm     atm bottom level pressure (pa)
c
            pgcm(i) = psurf(i)*sl1
c
c 2.3  wind     atm bottom level wind speed (m/s)
c      wind(i) = wind(i)
c
c 2.4  tgcm     atm bottom level temperature (kelvin)
c
            tgcm(i) = t1(i)
c
c 2.5  qgcm     atm bottom level specific humidity (kg/kg)
c
            qgcm(i) = q0(i)
c
c 2.6  zgcm     atm bottom level height above surface (m)
c
            zgcm(i) = -rd * tv1(i) * log(sl1) / g
c
c 2.7  flwds    downward longwave rad onto surface (w/m2)
c               (positive, not net longwave)
c
            flwds(i) = dlwflx(i)
c
c 2.8  sols     solar downward radiation (w m-2; positive, not net solar)
c
            sols(i) = dswflx(i)
         end if
      end do
c
c
c 3. soil parameters
c
c      binf ,ds   ,dsm  ,ws   ,cef  ,expt ,kst  ,dph  ,bub  ,qrt  ,
c      bkd  ,sld  ,wcr  ,wpw  ,silz ,snwz ,smr  ,smx  ,dphnd,smxnd,
c      expnd,bubnd,alpnd,betnd,gamnd
c
      do i = 1, im
         if(slimsk(i).eq.1.) then       
            do k = 1, msl
               expt(k,i)= expld(i,k)
               kst(k,i) = kstld(i,k)
               dph(k,i) = dphld(i,k)
               bub(k,i) = bubld(i,k)
               qrt(k,i) = qrtld(i,k)
               bkd(k,i) = bkdld(i,k)
               sld(k,i) = sldld(i,k)
               wcr(k,i) = wcrld(i,k)
               wpw(k,i) = wpwld(i,k)
               smr(k,i) = smrld(i,k)
               smx(k,i) = smxld(i,k)
            end do
            do k = 1, nsl
               dphnd(k,i) = dpnld(i,k)
               smxnd(k,i) = sxnld(i,k)
               expnd(k,i) = epnld(i,k)
               bubnd(k,i) = bbnld(i,k)
               alpnd(k,i) = apnld(i,k)
               betnd(k,i) = btnld(i,k)
               gamnd(k,i) = gmnld(i,k)
            end do
         endif
      enddo
c
c 4. vegetation parameters
c
c 4.1  nveg(imx2)              number of vegetation type in a grid cell
c 4.2  nvt(msub,imx2)          vegetation type
c 4.3  wt(msub,imx2)           fraction of grid cell covered by veg
c 4.4  rt(msl,msub,imx2)       fraction root in the soil layer
c 4.5  flai(msub,imx2)         leaf area index
c
      do i = 1, im
         if(slimsk(i).eq.1.) then       
            nveg(i) = 1

            if(sigmaf(i).gt.1.e-5.and.sigmaf(i).lt.(1.0-1.e-5))then
               if(vegtype(i).gt.0 .and. vegtype(i).lt.12) then
                  nveg(i) = 2
               end if
            end if
            if(nveg(i).eq.2) then
               nvt(1,i) = vegtype(i)
               nvt(2,i) = 12                ! vic bare soil
               wt(1,i)  = sigmaf(i)
               wt(2,i)  = 1.0 - sigmaf(i)
               do k = 1, msl
                  rt(k,1,i) = vrtld(i,k)
                  rt(k,2,i) = 0.
               end do
               flai(1,i) = flaild(i)
               flai(2,i) = 0.0
            else
               if(vegtype(i).lt.12) then
                  nvt(1,i) = vegtype(i)
               else
                  nvt(1,i) = 12
               endif
               wt(1,i)  = 1.0
               do k = 1, msl
                  rt(k,1,i) = vrtld(i,k)
               end do
               flai(1,i) = flaild(i)
            end if
         end if
      end do
c
c
c 5. modified variables 
c
c 5.1  cwt(md,msub,imx2)      canopy water (mm h2o)
c 5.2  csn(md,msub,imx2)      canopy snow  (mm h2o)
c 5.3  smc(msl,md,msub,imx2)  soil moisture content (mm)
c 5.4  sic(msl,md,msub,imx2)  soil ice content (mm)
c 5.5  tnd(nsl,md,msub,imx2)  temperature at soil nodes (k)
c 5.6  swq(md,msub,imx2)      snow water equivalent (mm h2o)
c 5.7  rsn(md,msub,imx2)      snow density (kg/m^3)
c 5.8  tsf(md,msub,imx2)      snow surface temperature (k)
c 5.9  tpk(md,msub,imx2)      snow pack temperature (k)
c 5.a  sfw(md,msub,imx2)      surface snow water equivalent (mm h2o)
c 5.b  pkw(md,msub,imx2)      snow pack snow water equivalent (mm h2o)
c 5.c  fmu(imx2)              precipitation fraction
c 5.d  lstsn(imx2)            number of time step since last snow
c
      do i = 1, im
         if(slimsk(i).eq.1.) then 
c            do jd = 1, md
c               cwt(jd,1,i) = canopy(i)/1000.0   ! mm -> m
c               csn(jd,1,i) = csno(i)            ! m
c            end do
            do nv = 1, nveg(i)
               do jd = 1, md
                  do m = 1, msl
                     smc(m,jd,nv,i) = smcld(i,m)*dph(m,i)*1000.0  ! mm
                     sic(m,jd,nv,i) = sicld(i,m)
                  end do
                  do n = 1, nsl
                     tnd(n,jd,nv,i) = stcld(i,n)
                  end do
                  swq(jd,nv,i) = sheleg(i)/1000.0    ! mm -> m
                  rsn(jd,nv,i) = rsno(i)

                  tsf(jd,nv,i) = tsfld(i)
                  tpk(jd,nv,i) = tpkld(i)
                  sfw(jd,nv,i) = sfwld(i)
                  pkw(jd,nv,i) = pkwld(i)
                  if(nv.eq.1)then
                     cwt(jd,nv,i) = canopy(i)/1000.0
                     csn(jd,nv,i) = csno(i)
                  else
                     cwt(jd,nv,i) = 0.0
                     csn(jd,nv,i) = 0.0
                  endif
               end do
            end do
         end if
         fmu(i) = 1.0
      end do
c
c
c    call vic lsm
c
#ifdef DBG
      print*,' --- in sfc_vic1 --- before sfuwvic'
      print*,'lat=',lat,' im =',im
#endif
c      print*,'lat=',lat,' im =',im
c
      do i = 1, im
c         if((lat.eq.1.and.(i.eq.113.or.i.eq.117)).or.
c     &      (lat.eq.3.and.(i.eq.85.or.i.eq.121)).or.
c     &      (lat.eq.6.and.i.eq.164).or.(lat.eq.9.and.i.eq.130))then

c         if(lat.eq.31.and.i.eq.98)then
c            print*,' --- in sfc_vic1 --- before sfuwvic'
c            print*,'lat=',lat,' im =',im,' i=',i

         if(slimsk(i) .eq. 1.) then
            sh(i) = 0.0
            lh(i) = 0.0
            gh(i) = 0.0
            ts(i) = 0.0
            albd(i) = 0.0
            ovflw(i) = 0.0
            bsflw(i) = 0.0
            sfall(i) = 0.0
c
#ifdef DBG
      print*,'before sfuwvic'
      print*, i, msl,          nsl,         msub,         md,
     1                 dtime,           mm,
     2                prc(i),      pgcm(i),      wind(i),    tgcm(i),
     3               qgcm(i),      zgcm(i),     flwds(i),    sols(i),
     4               binf(i),        ds(i),       dsm(i),      ws(i),
     5             expt(1,i),     kst(1,i),     dph(1,i),     cef(i),
     6              bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),
     7              wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),
     8              smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),
     9            expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),
     a            gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),
     b             rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),
     c          smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),
     d            rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),
     e            pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),
     f                 sh(i),        lh(i),        gh(i),      ts(i),
     g               albd(i),     ovflw(i),     bsflw(i),  snowmt(i),
     h             snowev(i)
#endif
c
            call sfuwvic(msl,          nsl,         msub,         md,
     1                 dtime,           mm,            i,        lat,
     2                prc(i),      pgcm(i),      wind(i),    tgcm(i),
     3               qgcm(i),      zgcm(i),     flwds(i),    sols(i),
     4               binf(i),        ds(i),       dsm(i),      ws(i),
     5             expt(1,i),     kst(1,i),     dph(1,i),     cef(i),
     6              bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),
     7              wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),
     8              smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),
     9            expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),
     a            gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),
     b             rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),
     c          smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),
     d            rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),
     e            pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),
     f                 sh(i),        lh(i),        gh(i),      ts(i),
     g               albd(i),     ovflw(i),     bsflw(i),  snowmt(i),
     h             snowev(i))
c
#ifdef DBG
            print *,'in sfc_vic1 i=',i,' lstsn=',lstsn(i)
            print *,'   tgcm=',tgcm(i),' ts=',ts(i)
            print *,'   tsf =',tsf(1,1,i),' tnd=',tnd(1,1,1,i)
c
      print*,'after sfuwvic'
      print*, i, msl,          nsl,         msub,         md,
     1                 dtime,           mm,
     2                prc(i),      pgcm(i),      wind(i),    tgcm(i),
     3               qgcm(i),      zgcm(i),     flwds(i),    sols(i),
     4               binf(i),        ds(i),       dsm(i),      ws(i),
     5             expt(1,i),     kst(1,i),     dph(1,i),     cef(i),
     6              bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),
     7              wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),
     8              smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),
     9            expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),
     a            gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),
     b             rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),
     c          smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),
     d            rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),
     e            pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),
     f                 sh(i),        lh(i),        gh(i),      ts(i),
     g               albd(i),     ovflw(i),     bsflw(i),  snowmt(i),
     h             snowev(i)
c
      print*,'- in sfc_vic1 - after sfuwvic ts (i) =',ts(i),sh(i)
#endif
         endif 

c         endif
      enddo
c
c 6.  update the modified variables
c
      do i = 1, im
         if(slimsk(i).eq.1.) then 
            canopy(i) = 0.
            csno(i)   = 0.
            do m = 1, msl
               smcld(i,m) = 0.0
               sicld(i,m) = 0.0
            end do
            do n = 1, nsl
               stcld(i,n) = 0.0
            end do
            sheleg(i)= 0.0
            rsno(i)  = 0.0
            tsfld(i) = 0.0
            tpkld(i) = 0.0
            sfwld(i) = 0.0
            pkwld(i) = 0.0

            do jd = 1, md           ! since just one-tile
               if(jd.eq.1) then
                  canopy(i) = canopy(i) + cwt(jd,1,i)*fmu(i)
                  csno(i)   = csno(i) + csn(jd,1,i)*fmu(i)
               else
                  canopy(i) = canopy(i) + cwt(jd,1,i)*(1-fmu(i))
                  csno(i)   = csno(i) + csn(jd,1,i)*(1-fmu(i))
               endif
            enddo

            do nv = 1, nveg(i)
               do jd = 1, md
                  if(jd.eq.1) then
                     wtld = fmu(i)*wt(nv,i)
                  else
                     wtld = (1.0-fmu(i))*wt(nv,i)
                  end if
c                  canopy(i) = canopy(i) + cwt(jd,nv,i)*wtld
c                  csno(i)   = csno(i) + csn(jd,nv,i)*wtld
c
                  do m = 1, msl
                     if(dph(m,i).gt.0) then
                        smc(m,jd,nv,i)=smc(m,jd,nv,i)/dph(m,i)/1000.0
                     else
                        smc(m,jd,nv,i) = 0.0
                     end if
                     smcld(i,m) = smcld(i,m) + smc(m,jd,nv,i)*wtld
                     sicld(i,m) = sicld(i,m) + sic(m,jd,nv,i)*wtld
                  end do
                  do n = 1, nsl
                     stcld(i,n) = stcld(i,n) + tnd(n,jd,nv,i)*wtld
                  end do

                  sheleg(i)= sheleg(i) + swq(jd,nv,i)*wtld
                  rsno(i)  = rsno(i) + rsn(jd,nv,i)*wtld
                  tsfld(i) = tsfld(i) + tsf(jd,nv,i)*wtld
                  tpkld(i) = tpkld(i) + tpk(jd,nv,i)*wtld
                  sfwld(i) = sfwld(i) + sfw(jd,nv,i)*wtld
                  pkwld(i) = pkwld(i) + pkw(jd,nv,i)*wtld
               end do
            end do
            sheleg(i) = sheleg(i) * 1000.0    ! m -> mm
            canopy(i) = canopy(i) * 1000.0    ! m -> mm
         end if
      end do
c
c   vic -> gsm/rsm: prepare variables for return to parent model
c
c 7. output (o):
c
c   return the following output fields to parent model
c
c  snowmt,snowev,gflux,
c  drain,evap,hflx,ep,
c  snowfl, runoff
c
      do i = 1, im
         gflux(i) = 0.0              ! initialize nonland grid cell
         snowfl(i) = 0.0
         runoff(i) = 0.0
         drain(i) = 0.0

         if(slimsk(i) .eq. 1.) then
c
c   tskin  land surface temperature (k)
c
            tskin(i) = ts(i)
c
c   evap   actual latent heat flux (w m-2: positive, if upward from sfc)
c
            evap(i) = - lh(i)
c
c   hflx       sensible heat flux (w m-2: positive, if upward from sfc)
c
            hflx(i) = - sh(i)
c
c   gflux      soil heat flux (w m-2: negative if downward from surface)
c
            gflux(i) = gh(i)
c
c   ep         potential evaporation (w m-2)
c
* -- temporary, use osulsm2 to compute ep ------------------------------
c
c  compute potential evaporation for land
c
            t2 = ts(i) * ts(i)
            t4 = t2 * t2
c
c  rcap = fnet - sigma t**4 + gflx - rho cp ch v (t1-theta1)
c
            rcap = (1.0-albd(i))*sols(i) - sigma * t4 + gflux(i)
     &             -rch(i) * (ts(i) - theta1(i))
c
c  rsmall = 4 sigma t**3 / rch + 1
c
            rsmall = 4.*sigma*ts(i)*t2/rch(i) + 1.
c
c  delta = l / cp * dqs/dt
c
            delta = elocp * eps * hvap * qs1(i) / (rd * t2)
c
            ep(i) = elocp*rsmall*rch(i)*(qs1(i)-q0(i))+rcap*delta
            ep(i) = ep(i) / (rsmall + delta)
c
* -- end of computing ep -----------------------------------------------
c
c   runoff     surface runoff (mm s-1), not infiltrating the surface
c   drain    subsurface runoff (mm s-1), drainage out bottom
c
            snowfl(i) = sfall(i)                 ! m/time step
            runoff(i) = ovflw(i)/dtime
            drain(i)  = bsflw(i)/dtime
         endif
      enddo
c
c   dm       ratio of actual/potential evap (dimensionless)
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
