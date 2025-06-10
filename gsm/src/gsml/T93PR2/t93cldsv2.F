!----------------------------------------------------------------------
! original file is t93clds.F90 from ECPC SCM by S.Iacobellis
!
!   may. 2006     modified by A.Shimpo
!     * convert freeform format to fixed format
!     * modified subroutine name "t93clds" to "t93cldsv2"
!     * setting "define TUNE1" : same as aug.2005 modification in t93clds.F
!       ** add calculation of tropopause level = itrop
!       ** modify rhcrit when sgk(i)<=0.65.
!          rhcrit=0.80 -> 0.80+0.20*(0.65-sgk(i))/(0.65-sgk(itrop))
!       ** dqsdt=dqsdpma*(wpa(i)+0.05*g*cmf(i)) + ...
!         ->dqsdt=dqsdpma*(wpa(i)+g*cmf(i)) + ...
!       ** modify e(i) when lw(i)-e(i)*delt < 0
!       ** modify dcoef
!----------------------------------------------------------------------
!**********************************************************************
!!      subroutine t93clds(sgk   ,p     ,pk    ,t     ,th    ,
      subroutine t93cldsv2(sgk   ,p     ,pk    ,t     ,th    ,
     1                   q     ,clw   ,ciw   ,rho   ,dlf   ,du    ,
     2                   z     ,dtsur ,dqsur ,solin ,wpa   ,
     3                   cmf   ,hrtr  ,hrsl  ,dthor ,hpbl  ,
     4                   fqs   ,fss   ,pref  ,delt  ,a     ,
     5                   dtliw ,dqliw ,prec  ,precl ,precf ,iflag,
     6                   sacv  ,sac   ,sabl  ,disp  ,        
     7                                                       lfrc )
!!     6                   sacv  ,sac   ,sabl  ,disp  ,paut  ,
!!     7                   pcll  ,pacc  ,peva  ,pflx  ,psub  , lfrc )

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use ppgrid,        only: nl => pver, nl1 => pverp
      use physconst,     only: l => latvap,  g => gravit, cpa => cpair,
     1                         gam => cappa, ra => rair,  lf => latice
      implicit none

!!#define TUNE1  ! tuning in aug, 2005
 
!	cloud and liquid water parameterization of tiedtke, 1993
!
!	written to be used in CAM3 and is called once for each column.
!	in this subroutine, layer number increases from the ground up, so
!	cam3 arrays must be inverted prior to call. 

!	updated Dec, 2005 to include new precipitation physics (s.r. precip2)

!	input variables

      real(r8), intent(in) :: delt		! timestep (seconds)
      real(r8), intent(in) :: sgk(nl)		! sigma coordinates
      real(r8), intent(in) :: p(nl1)		! interface pressure (Pa)
      real(r8), intent(in) :: pk(nl)		! midlayer pressure (Pa)
      real(r8), intent(in) :: t(nl)		! temperature (deg K)
      real(r8), intent(in) :: th(nl)		! potential temperature (deg K)
      real(r8), intent(in) :: q(nl)		! specific humidity (kg/kg)
      real(r8), intent(in) :: rho(nl)		! density (kg/m**3)
      real(r8), intent(in) :: dlf(nl)		! detrained cloud water/ice from convection (kg/kg/sec)
      real(r8), intent(in) :: du(nl)		! mass detrainment rate (1/sec)
      real(r8), intent(in) :: z(nl1)		! interface height (meters)
      real(r8), intent(in) :: dtsur(nl)		! temperature change due to diffusion and surface flux (degK/sec)
      real(r8), intent(in) :: dqsur(nl)		! moisture change due to diffusion and surface flux (kg/kg/sec)
      real(r8), intent(in) :: solin(nl)		! net solar absorbtion (W/m**2)
      real(r8), intent(in) :: wpa(nl)		! large-scale vertical velocity (Pa/sec)
      real(r8), intent(in) :: cmf(nl)		! cumulus mass flux (kg/m**2/sec)
      real(r8), intent(in) :: hrtr(nl)		! ir heating rate (degK/sec in CAM3 application)
      real(r8), intent(in) :: hrsl(nl)		! solar heating rate (degK/sec in CAM3 application)
      real(r8), intent(in) :: dthor(nl)		! horiz. adv. heating (degK/timestep) not used in CAM3 application
      real(r8), intent(in) :: hpbl		! boundary layer height
      real(r8), intent(in) :: fqs		! surface latent heat flux
      real(r8), intent(in) :: fss		! surface sensisble heat flux
      real(r8), intent(in) :: pref		! reference pressure (mb)
      real(r8), intent(in) :: lfrc		! land fraction

!	input/output variables

      real(r8), intent(inout) :: a(nl)		! cloud fraction
      real(r8), intent(inout) :: clw(nl)	! cloud liquid water (kg/kg) (grid average)
      real(r8), intent(inout) :: ciw(nl)	! cloud ice water (kg/kg) (grid average)
 
!	output variables

      real(r8), intent(out) :: prec		! layer (stratiform) precipitation
      real(r8), intent(out) :: precl		! liquid layer (stratiform) precipitation
      real(r8), intent(out) :: precf		! frozen layer (stratiform) precipitation
      real(r8), intent(out) :: dtliw(nl)	! temperature change (deg K)
      real(r8), intent(out) :: dqliw(nl)	! humidity change (kg/kg)

!	locally computed variables

      real(r8) dp(nl)		! layer thickness (Pa)
      real(r8) dz(nl)		! layer thickness (meters)
      real(r8) lw(nl)		! total cloud water (liquid+ice) kg/kg (grid average)
      real(r8) scv(nl)		! cloud water source due to convection (kg/kg/sec)
      real(r8) scv_l, scv_i
      real(r8) sbl(nl)		! cloud water source due to stratocumulus (boundary layer) clouds
      real(r8) sbl_l, sbl_i
      real(r8) sacv(nl)		! cloud source due to convection (fraction/sec)
      real(r8) sac(nl)		! cloud source due to stratiform cloud (fraction/sec)
      real(r8) sabl(nl)		! cloud source due to stratocumulus (boundary layer) clouds (fraction/sec)
      real(r8) cond(nl)		! condensation rate in stratiform clouds
      real(r8) cond_l, cond_i
      real(r8) e(nl)		! evaporation rate of cloud water
      real(r8) e_l, e_i
      real(r8) disp(nl)		! cloud disipation rate (fraction/sec)
      real(r8) paut(nl)		! conversion of cloud water to precip due to autoconversion (g/m**2/sec)
      real(r8) pcll(nl)		! conversion of cloud water to precip due to collection (g/m**2/sec)
      real(r8) pacc(nl)		! conversion of cloud water to precip due to accretion (g/m**2/sec)
      real(r8) peva(nl)		! evaporation of liquid precipitation (g/m**2/sec)
      real(r8) pflx(nl)		! flux divergence of ice particles (g/m**2/sec)
      real(r8) psub(nl)		! sublimation of frozen precipitation (g/m**2/sec)
      real(r8) fs(nl1)		! turbulent heat flux at interfaces (W/m**2)
      real(r8) wpthp(nl1)	! turbulent heat flux at interfaces (degK m/sec)
      real(r8) fq(nl1)		! turbulent moisture flux at interfaces (W/m**2)
      real(r8) wpqp(nl1)	! turbulent moisture flux at interfaces (degK m/sec)
      real(r8) dq_ent(nl)	! rate of humidity change due to entrainment of stratocumulus 
      real(r8) dt_ent(nl)	! rate of temperature change due to entrainment of stratocumulus 
      real(r8) dl_ent(nl)	! rate of liquid water change due to entrainment of stratocumulus 
      real(r8) dlent_l, dlent_i
      real(r8) dtmelt(nl)	! temperature change due to melting of frozen precipitation (deg K)
      real(r8) dtprec(nl)	! temperature change due to melting/freezing of precipitation (deg K)
      real(r8) fliqex(nl)	! existing liquid water fraction
      real(r8) ficeex(nl)	! existing ice water fraction
      real(r8) dtdt_p(nl)	! change in temperature due to precip processes (deg K/sec)
      real(r8) dqdt_p(nl)	! change in specific humidity due to precip processes (kg/kg/sec)
      real(r8) dqldt_p(nl)	! change in cloud liquid water due to precip processes (kg/kg/sec)
      real(r8) dqidt_p(nl)	! change in cloud ice water due to precip processes (kg/kg/sec)
      real(r8) fliq		! real(r8) :: function
      real(r8) fl		! liquid fraction of cloud water (from prescribed function fliq)
      real(r8) sablmx		! maximum allowed value of sabl
      real(r8) plcl		! pressure at LCL
      real(r8) tlcl		! temperature at LCL
      real(r8) the2		! equivalent potential temperature of layer 2 air
      real(r8) thes		! saturated equivalent potential temperature of a given layer  
      real(r8) fqbase		! turbulent moisture transport at base of boundary layer cloud
      real(r8) q0bl		! humidity at base of boundary layer cloud
      real(r8) qtbl		! humidity at top of boundary layer cloud
      real(r8) qstbl		! sat. humidity at top of boundary layer cloud
      real(r8) ttbl		! temperature at top of boundary layer cloud
      real(r8) lwtop		! cloud water/ice near top of b.l. cloud
      real(r8) prate		! precipitation rate (kg/m**2/sec) (stratiform)
      real(r8) lprate		! liquid precipitation rate (kg/m**2/sec) (stratiform)
      real(r8) fprate		! frozen precipitation rate (kg/m**2/sec) (stratiform)
      real(r8) lwubl		! liquid water in b.l. updraft (kg/kg)
      real(r8) lwdbl		! liquid water in b.l. downdraft (kg/kg)
      real(r8) thsf		! potential temperature near surface
      real(r8) ptop		! pressure near top of b.l. cloud
      real(r8) pmidc		! pressure near midpoint of b.l. cloud
      real(r8) pbase		! pressure near base of b.l. cloud
      real(r8) r1, r2		! interpolation dummy variables
      real(r8) ratio
      real(r8) rhows
      real(r8) tflxnet
      real(r8) qflxnet
      real(r8) denlev
      real(r8) zoverh
      real(r8) pmidmx
      real(r8) thovert
      real(r8) wpthep0
      real(r8) wpthept
      real(r8) sumabs
      real(r8) sumrho
      real(r8) rhomx
      real(r8) rshmrso
      real(r8) rscmrso
      real(r8) de1, de2, de3, de4
      real(r8) wpqwp0
      real(r8) wpqwpt
      real(r8) wpthvpb
      real(r8) wpthvpa
      real(r8) wpthvpm
      real(r8) zmida
      real(r8) zmidb
      real(r8) del_sv
      real(r8) del_s
      real(r8) del_q
      real(r8) del_l
      real(r8) we
      real(r8) wpsph
      real(r8) wpqph
      real(r8) wplph
      real(r8) ds_entb
      real(r8) ds_enta
      real(r8) dqsdpma			! change of qs along moist adiabat
      real(r8) dqsdtp			! change of qs w.r.t. temperature
      real(r8) th1, th2			! used to calculate dqsdpma and dqsdtp
      real(r8) p1, p2, p3		! used to calculate dqsdpma and dqsdtp
      real(r8) t1, t2, t3		! used to calculate dqsdpma and dqsdtp
      real(r8) qs1, qs2, qs3		! used to calculate dqsdpma and dqsdtp
      real(r8) dtdtd			! temperature tendency due to diabatic processes
      real(r8) dqsdt			! saturation humidity tendency
      real(r8) c1			! condensation in existing clouds
      real(r8) c2			! condensation in new clouds
      real(r8) rhcrit			! critical RH for condensation to occur
      real(r8) qs
      real(r8) rhlay
      real(r8) sacmx
      real(r8) e1			! evaporation due to descent
      real(r8) e2			! evaporation due to turbulent mixing
      real(r8) pratenew			! new precipitation rate calculated if analytic integration is used
      real(r8) lpratenew		! new liquid precipitation rate calculated if analytic integration is used
      real(r8) fpratenew		! new frozen precipitation rate calculated if analytic integration is used
      real(r8) dqfix			! correction term
      real(r8) aa, bb			! terms used to calculate new cloud amount (analytic integration)
      real(r8) cc, cc_l, cc_i		! terms used to calculate new cloud water (analytic integration) 
      real(r8) dd, dd_l, dd_i		! terms used to calculate new cloud water (analytic integration) 
      real(r8) oldlw			! used to adjust precipitation if analytic integration is used
      real(r8) gpnew
      real(r8) tw_start, tw_end		! total water amounts at start and end of routine (used to balance water budget)
      real(r8) te_start, te_end		! total energy amounts at start and end of routine (used to balance water budget)
      real(r8) dtcol,dqcol, wvcol, wlcol, wicol, secol
      real(r8) dtcolabs, dqcolabs, tcol, qcol
      real(r8) dtadjr,dqadjr,dqbef,dlfcol
      real(r8) clwb(nl), ciwb(nl)
      integer lbase		! model layer containing the base of any b.l. cloud
      integer layb
      integer layt
      integer nlay
      integer i, k
      integer lay, lev
      integer istart
      integer ifd		! flag determines if analytic integration is used

      real(r8) aa1(nl)		! term used for analytic integration
      real(r8) aa2(nl)		! term used for analytic integration
      real(r8) aa3(nl)		! term used for analytic integration
      real(r8) bb1(nl)		! term used for analytic integration
      real(r8) dd1(nl)		! term used for analytic integration
      real(r8) :: qsat		! real(r8) :: function

      real(r8) :: sumtot,sumfix,summlt,sumcnd,sume,sumep,sumepr,sumgpr
      real(r8) :: sumqv0,sumql0,sumqi0
      real(r8) :: sumqv1,sumql1,sumqi1
      real(r8) :: sumqtot,sumtest
      integer :: iflag

      real(r8) :: sumq,sumql,sumqi

#ifdef TUNE1
      integer :: itrop
      real(r8) :: dthdp
#endif

#ifndef TUNE1
      real, parameter :: dcoef = 1.0e-06	! diffusion coefficient for evaporation of cloud water from turbulent mixing (1/sec)
#else
      real, parameter :: dcoef = 2.0e-06	! diffusion coefficient for evaporation of cloud water from turbulent mixing (1/sec)
#endif

!	initialize necessary variables

      do 1 i=1,nl
      dq_ent(i)=0.0
      dt_ent(i)=0.0
      dl_ent(i)=0.0
      sabl(i)=0.0
      sbl(i)=0.0
      fq(i)=0.0
      fs(i)=0.0
      wpthp(i)=0.0
      wpqp(i)=0.0
      dtmelt(i)=0.0
      dtprec(i)=0.0
      aa1(i)=0.0
      aa2(i)=0.0
      aa3(i)=0.0
      bb1(i)=0.0
      dd1(i)=0.0
 1    continue

      do 5 i=1,nl
      if (clw(i)+ciw(i) .eq. 0.0) then
        fliqex(i) = 0.0
        ficeex(i) = 0.0
      else
        fliqex(i) = clw(i) / (clw(i) + ciw(i))
        ficeex(i) = ciw(i) / (clw(i) + ciw(i))
      endif
      lw(i)=clw(i)+ciw(i)
      dp(i)=p(i)-p(i+1)
      dz(i)=z(i)-z(i+1)			!caution: negative values here
 5    continue

#ifdef TUNE1
!
!   tropopause level   ( see cldjmsn.F )
!     i<=itrop : troposphere, i>itrop : stratosphere
!
      itrop=nl
      do i=2,nl
        if (pk(i).le.450.e2) then
          dthdp=(th(i+1)-th(i))/(pk(i+1)-pk(i))
          if (dthdp.lt.-0.25e-2) then
            itrop=i-1
            exit
          endif
        endif
      end do
#endif

      secol=0.0
      wvcol=0.0
      wlcol=0.0
      wicol=0.0
      dlfcol=0.0
      do 6 i=1,nl
        secol = secol + cpa*t(i)*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
        dlfcol = dlfcol + dlf(i)*delt*dp(i)/g
 6    continue
      te_start = secol + (l+lf)*wvcol + lf*wlcol + lf*dlfcol
      tw_start = wvcol + wlcol + wicol + dlfcol
!     if(iflag.eq.1) then
!       write(6,*)'*****************************************'
!       write(6,*)'INSIDE T93CLDS, AT START'
!       write(6,*)'wvcol = ', wvcol
!       write(6,*)'wlcol = ', wlcol
!       write(6,*)'wicol = ', wicol
!       write(6,*)'dlfcol= ', dlfcol
!       write(6,*)'secol = ', secol
!       write(6,*)'tw_start=',tw_start
!       write(6,*)'te_start=',te_start
!     endif 

!	1.0) cloud water and cloud fraction terms from cumulus convection

      do 10 i=1,nl
      scv(i)=dlf(i)				! rate of change in liquid water - convection (kg/kg/sec) (dlf includes both Zhang and Hack)
      sacv(i)=(1.0-a(i))*du(i)			! rate of change in cloud cover - convection (du units = 1/sec)
      aa1(i)=du(i)				! save this term for analytic integration
 10   continue

!	2.0) cloud water and cloud fraction terms from
!	     boundary layer clouds (stratocumulus)

!	2.1) find the layer which contains the top of the 
!	     planetary boundary layer. lbase may be reset below.

      do 20 i=1,nl
      if(hpbl.ge.z(i).and.hpbl.lt.z(i+1)) then
        lbase=i
        goto 21
      endif
 20   continue
 21   lbase=max(lbase,2)			

!	2.2) use definitions in Betts and Albrecht (jas, 1987, p83-99) 
!	     to establish the cloud region within the convective 
!	     boundary layer (cbl). the cloud region is defined as those
!	     layers where the saturated equivalent potential
!	     temperature (thes) is less than the equivalent potential
!            temperature (the) in the layer where the lifting process 
!            starts (here defined as layer 2).

      call lcl(pk  ,th  ,q   ,2   ,plcl ,pref )		!find lcl for using layer 2 as starting level
      tlcl=th(2)*((plcl/pref)**gam)			!temperature at lcl
      the2=th(2)*exp(l*q(2)/(cpa*tlcl))			!equivalent potential temp of layer 2 air

      istart=0
      layb=0
      layt=0
      do 3838 k=2,nl

      if(sgk(k).lt.0.6) then				!if reach this high before seeing inversion
        nlay=10						!then set nlay=10 and exit this loop because
        goto 3839					!either the CBL is non-existant or the cloud
      endif						!layer extends across more than one model layer

      thes=th(k)*exp(l*qsat(t(k),pk(k))/(cpa*t(k)))	!saturation equivalent potential temperature of layer k

      if(istart.eq.0.and.thes.le.the2) then		!enter this block if we are encountering positive
        istart=1					!buoyancy (thes<the2) for the first time.
        layb=k						!set layb as the bottom of this cloud
      endif

      if(istart.eq.1.and.thes.gt.the2) then		!enter this block if we have reached the top of
        nlay=k-layb					!the cloud layer defined as (thes>the2). set nlay 
        goto 3839					!to the number of model layers within this cloud
      endif						!and exit loop 3838

 3838 continue
 3839 if(nlay.gt.1) goto 401				!either no condensation or cloud spans more than one model layer
      lbase=layb
      if(sacv(lbase).gt.0.0.or.scv(lbase).gt.0.0) goto 401	!skip section if convection present 
      fqbase=dqsur(lbase)*dp(lbase)/g				!approximate moisture flux thru base of cloud

!	2.3) calculate cloud mass flux at base of cloud (rhows).
!	     use interpolation to obtain cloud top values.

      q0bl=q(1)					!humidity near surface
      r1=(pk(lbase)-p(lbase+1))/(pk(lbase)-pk(lbase+1))
      r2=1.0-r1
      qtbl=r1*q(lbase+1)+r2*q(lbase)		!humidity near cloud top
      ttbl=r1*t(lbase+1)+r2*t(lbase)		!temperature near cloud top
      qstbl=qsat(ttbl,p(lbase+1))		!sat. humidity near top of b.l. cloud

      lwtop=lw(lbase+1)				!cloud water near top of b.l. cloud (use grid-mean value)
      ratio=q0bl/(a(lbase)*(qstbl+lwtop)+qtbl*(1-a(lbase)))
      if(ratio.lt.1.00001) then
        rhows=0.0
      else
        rhows=fqbase/(q0bl-(a(lbase)*(qstbl+lwtop)+qtbl*(1-a(lbase))))
      endif
      rhows=max(0.0,rhows)

!     if(iflag.eq.1) then
!       write(6,*)'lbase=',lbase
!       write(6,*)'dp(lbase)=',dp(lbase)
!       write(6,*)'dqsur(lbase)=',dqsur(lbase)
!       write(6,*)'fqbase=',fqbase
!       write(6,*)'rhows=',rhows
!       write(6,*)'q0bl=',q0bl
!       write(6,*)'a(lbase)=',a(lbase)
!       write(6,*)'qstbl=',qstbl
!       write(6,*)'lwtop=',lwtop
!       write(6,*)'qtbl=',qtbl
!       write(6,*)'ratio=',ratio
!     endif

!	2.4) determine liquid water in updrafts (lwubl) and 
!	     downdrafts (lwdbl). evaluate these at midpoint of
!	     lbase

!	updraft properties

      thsf=th(1)				!potential temp. near surface

!	properties near top, middle and base of b.l. cloud

      lwtop=lw(lbase)
      ptop=p(lbase+1)
      pmidc=pk(lbase)
      pbase=p(lbase)

      call getchars(thsf  ,ttbl  ,pbase ,pmidc ,ptop  ,
     1              lwubl ,lwdbl ,l     ,cpa   ,gam   ,
     2              pref  )

      lwdbl=max(lwtop+lwdbl,0.0)

!	2.5) now calculate cloud water and cloud fraction rates.
!	     impose two constraints. first, limit sbl*delt to 10% of
!	     available q. Second limit sabl not to exceed (1-a)/delt or
!	     2.78e-04 sec-1 (clear to overcast in less than 1hr).
!	     if exceeded, set sabl=(1-a)/delt and recalculate rhows
!	     and sbl.

      sbl(lbase)=-rhows*(lwubl-a(lbase)*lwdbl)/(rho(lbase)*dz(lbase))
      if(sbl(lbase)*delt .gt. 0.10*q(lbase)) then
        sbl(lbase) = 0.10*q(lbase)/delt
        rhows = sbl(lbase)/(lwubl-a(lbase)*lwdbl)/(rho(lbase)*dz(lbase))
      endif

      sabl(lbase)=-rhows*(1.-a(lbase))/(rho(lbase)*dz(lbase))
      sablmx=min( (1.0-a(lbase))/delt , 2.78e-04 )
      if(sabl(lbase).gt.sablmx) then
        sabl(lbase)=sablmx
        rhows=-sabl(lbase)*rho(lbase)*dz(lbase)/(1.-a(lbase))
        sbl(lbase)=-rhows*(lwubl-a(lbase)*lwdbl)/(rho(lbase)*dz(lbase))
      endif
      aa2(lbase)=-rhows/(rho(lbase)*dz(lbase))	!save for analytic integration

!	3.0) entrainment terms at top of stratocumulus cloud

!	3.1) determine turbulent fluxes of theta_e at surface (level=1) 
!            and top of cloud layer (level=lbase+1)

      fs(1)=fss
      fq(1)=fqs
      do 30 lay=1,lbase
      tflxnet=dtsur(lay)*cpa*dp(lay)/g		!net flux needed for layer
      qflxnet=dqsur(lay)*l*dp(lay)/g		!net flux needed for layer
      fs(lay+1)=fs(lay)-tflxnet
      fq(lay+1)=fq(lay)-qflxnet
 30   continue

      do 31 lev=1,lbase+1
      if(lev.eq.1) denlev=rho(1)
      if(lev.ne.1) denlev=0.5*(rho(lev-1)+rho(lev))
      wpthp(lev)=fs(lev)/(denlev*cpa)		!turbulent flux of th at level
      wpqp(lev)=fq(lev)/(denlev*l)		!turbulent flux of q at level
 31   continue
      zoverh=z(lbase)/z(lbase+1)				!here assume lcl is base of cloud layer
      pmidmx=0.5*(p(lbase+1)+p(1))				!pressure at middle of mixed layer
      thovert=(pref/pmidmx)**gam				!theta/temp for mixed layer
      wpthep0=wpthp(1)+thovert*l*wpqp(1)/cpa			!turbulent flux of theta_e at surface
      wpthept=wpthp(lbase+1)+thovert*l*wpqp(lbase+1)/cpa	!turbulent flux of theta_e at top of b.l. cloud layer

!	3.2) now calculate the solar flux terms from deardorff, 1976

      sumabs=0.0
      sumrho=0.0
      do 32 lay=1,lbase
      sumabs=sumabs+solin(lay)				!sum solar absorption in mixed layer (W/m**2)
      sumrho=sumrho+rho(lay)*dp(lay)
 32   continue
      rhomx=sumrho/(p(1)-p(lbase+1))			!mean density in mixed layer
      rshmrso=-sumabs/(rhomx*cpa)
      rscmrso=-(sumabs-solin(lbase))/(rhomx*cpa)

!	3.3) calculate the turbulent buoyancy flux below (wpthvpb) and
!            above (wpthvpa) cloud base. then get mean buoyancy flux
!            in mixed layer (wpthvpm)

      de1=1.005							!values from Deardorff, 1976
      de2=8.3							!value from Deardorff, 1976
      de3=0.05							!value from Deardorff, 1976
      de4=0.5							!value from Deardorff, 1976
      wpqwp0=0.0						!assume no flux of cloud water
      wpqwpt=0.0						!assume no flux of cloud water
      wpthvpb=de1*((1.-zoverh)*wpthep0 + zoverh*wpthept + 
     1            zoverh*thovert*rshmrso - thovert*rscmrso) -
     2            de2*thsf*((1.-zoverh)*wpqwp0 + zoverh*wpqwpt)

      wpthvpa=de4*((1.-zoverh)*wpthep0 + zoverh*wpthept +
     1           zoverh*thovert*rshmrso - thovert*rscmrso) -
     2           thsf*((1.-zoverh)*wpqwp0 + zoverh*wpqwpt) 

      wpthvpm=(p(1)-p(lbase))*wpthvpb + (p(lbase)-p(lbase+1))*wpthvpa
      wpthvpm=wpthvpm/(p(1)-p(lbase+1))


!	3.4) calculate entrainment velocity (we) and changes across
!	     inversion of the virtual dry static energy (del_sv), the
!	     dry static energy (del_s), the humidity (del_q), and the
!	     liquid water (del_l)

      zmida=0.5*(z(lbase+1)+z(lbase+2))
      zmidb=0.5*(z(lbase)+z(lbase+1))
      del_sv=cpa*t(lbase+1)*(1.+0.608*q(lbase+1))+g*zmida - 
     1       cpa*t(lbase)*(1.+0.608*q(lbase))+g*zmidb
      del_s=cpa*t(lbase+1)+g*zmida - cpa*t(lbase)+g*zmidb
      del_q=q(lbase+1) - q(lbase)
      del_l=lw(lbase+1) - lw(lbase)
      we=0.5*wpthvpm/del_sv

!	3.5) calculate the changes in temperature, humidity, and 
!	     liquid water due to entrainment.

      wpsph=-we*del_s
      wpqph=-we*del_q
      wplph=-we*del_l

!	note: switch sign since dz's are negative

      ds_entb=wpsph/dz(lbase)
      ds_enta=-wpsph/dz(lbase+1)
!     dt_ent(lbase)=ds_entb/cpa			!units=(deg K/sec)
!     dt_ent(lbase+1)=ds_enta/cpa

!     dq_ent(lbase)=wpqph/dz(lbase)		!units=(kg/kg/sec)
!     dq_ent(lbase+1)=-wpqph/dz(lbase+1)

!     dl_ent(lbase)=wplph/dz(lbase)		!units=(kg/kg/sec)
!     dl_ent(lbase+1)=-wplph/dz(lbase+1)

!	use this formulation for changes in layer lbase+1 to ensure energy/water balance

!     dt_ent(lbase+1)=-dt_ent(lbase)*dp(lbase)/dp(lbase+1)
!     dq_ent(lbase+1)=-dq_ent(lbase)*dp(lbase)/dp(lbase+1)
!     dl_ent(lbase+1)=-dl_ent(lbase)*dp(lbase)/dp(lbase+1)

!	4.0) formation of stratiform clouds

!	4.1) calculate the change of qs along moist adiabat (dqsdpma)
!  	     and change of qs with respect to Temp (dqsdtp)

 401  do 40 i=1,nl

      p1=pk(i)
      p2=pk(i)+1000.0
      t1=t(i)
      th1=t1*(pref/p1)**gam
      call getthc2s(p1   ,p2   ,th1  ,th2  ,l    ,
     1              cpa  ,gam  ,pref )
      t2=th2*(p2/pref)**gam
      qs1=qsat(t1,p1)
      qs2=qsat(t2,p2)
      dqsdpma=(qs2-qs1)/(p2-p1)
      p3=pk(i)
      t3=t(i)+0.5
      qs3=qsat(t3,p3)
      dqsdtp=(qs3-qs1)/(t3-t1)

!	4.2) calculate temperature tendency due to diabatic 
!	     processes (dtdtd)

      dtdtd = hrtr(i) + hrsl(i) + dtsur(i) + dt_ent(i) 

!	4.3) calculate saturation humidity tendency (dqsdt)

#ifndef TUNE1
      dqsdt=dqsdpma*(wpa(i)+0.05*g*cmf(i)) + dqsdtp*dtdtd       !assume 0.05 ratio between cmf and
#else
      dqsdt=dqsdpma*(wpa(i)+     g*cmf(i)) + dqsdtp*dtdtd       !assume 0.05 ratio between cmf and
#endif

!********************************************************
!	new handling of dqsdt when negative (may 2004, v2.75)

      if(dqsdt.lt.0.0) then
        dqsdt = dqsdt/(1.0+a(i)*dqsdtp*l/cpa)
      endif

!********************************************************
!	4.4) determine condensation in existing clouds (c1)

      if(dqsdt.lt.0.0) then
        c1=-a(i)*dqsdt						!explicit 
      else
        c1=0.0
      endif

!	4.5) calculate the formation rate of new stratiform
!	     clouds (sac). make sure sac does not exceed 
!	     limit of (1 - a)/delt. also check if q(i) > qs.
!	     if so, set sac to (1 - a)/delt.

      if(sgk(i).gt.0.65) then
        rhcrit=0.80+0.20*(sgk(i)-0.65)/(sgk(lbase)-0.65)
      else
#ifndef TUNE1
!       rhcrit=0.80+0.20*(0.65-sgk(i))/(0.65-0.10)		!sam's experiment 
        rhcrit=0.80
#else
        rhcrit=0.80+0.20*(0.65-sgk(i))/(0.65-sgk(itrop))
#endif
      endif
      rhcrit=max(rhcrit,0.0)
      rhcrit=min(rhcrit,1.0)
!      rhcrit=min(rhcrit,0.95)

      qs=qsat(t(i),pk(i))
      rhlay=q(i)/qs
      if(dqsdt.lt.0.0.and.rhlay.ge.rhcrit.and.qs.gt.q(i)) then
        sac(i)=-(1.0-a(i))*dqsdt/(qs-q(i))
        aa3(i)=-dqsdt/(qs-q(i))			!save for analytic integration
      else
        sac(i)=0.0
        aa3(i)=0.0				!save for analytic integration
      endif
      sacmx=(1.0-a(i))/delt
#ifndef TUNE1
      if((sac(i).gt.sacmx).or.(q(i).ge.qs)) then
#else
      if((sac(i).gt.sacmx).or.
     1   ( (dqsdt.lt.0.0).and.(q(i).ge.qs) ) ) then
#endif
        sac(i)=sacmx
        aa3(i)=1.0/delt				!save for analytic integration
      endif
#ifdef TUNE1
!
!   special treatment : dqsdt=0 above tropapause (option)
!
      if (dqsdt.lt.0.0.and.i.gt.itrop) then   ! above tropapause
        dqsdt=0.0
        c1=0.0
        sac(i)=0.0
        aa3(i)=0.0
      endif
#endif

!	4.6) calculate generation of cloud water from new
!	     clouds (c2). add to c1 to get total condensation (cond)

      if(dqsdt.lt.0.0) then
        c2=-delt*sac(i)*dqsdt			!explicit
      else
        c2=0.0
      endif

!     if(iflag.eq.1.and.i.eq.12) write(6,*)'i,c1,c2=',i,c1,c2
      cond(i)=c1+c2

      if(delt*cond(i) .gt. 0.25*q(i)) then		! limit amount of condensation to 25% of available q in any timestep
        cond(i) = 0.25*q(i)/delt
      endif

!     if(q(i)-delt*cond(i) .lt. 1.0e-10) then
!       cond(i)=(q(i)-1.0e-10)/delt
!       if(iflag.eq.1.and.i.eq.24) write(6,*)'after check, c1,c2,delt=',c1,c2,delt
!     endif

!	5.0) evaporation and dissipation of clouds

!	5.1) evaporation due to large-scale descent and 
!	     diabatic heating

      if(dqsdt.gt.0.0) then
        e1=a(i)*dqsdt						!explicit
      else
        e1=0.0
      endif

!	5.2) evaporation due to turbulent mixing

      e2=a(i)*dcoef*(qs-q(i))
      e2=max(e2,0.0)
      e(i)=e1+e2

!     if(iflag.eq.1.and.i.le.2) then
!       write(6,*)'before check, i=',i
!       write(6,*)'e1, e2, e=',e1,e2,e(i)
!     endif
      if(lw(i)-e(i)*delt .lt. 0.0) then
#ifndef TUNE1
        e1 = e1 * lw(i)/delt/e(i)
        e2 = e2 * lw(i)/delt/e(i)
        e(i) = lw(i)/delt
#else
        e(i) = lw(i)/delt
        e2 = min(e2, e(i))
        e1 = e(i) - e2
#endif
      endif
!     if(iflag.eq.1.and.i.le.2) then
!       write(6,*)'after check, i=',i
!       write(6,*)'e1, e2, e=',e1,e2,e(i)
!     endif

!	5.3) cloud dissipation. 
!		if lw > 0 set dissipation rate (eq 32, tiedtke 93)
!		if lw = 0 set dissipation to 0. if lw=0 at end of
!		          routine, a will be set to 0

      if(lw(i).gt.0.0) then		
        disp(i)=e2*a(i)/lw(i)
        bb1(i)=e2/lw(i)			!save for analytic integration
      else				
        disp(i)=0.0
        bb1(i)=0.0			!save for analytic integration
      endif

 40   continue

      sumtest=0.0
      sumtot=0.0
      sume=0.0
      sumep=0.0
      summlt=0.0
      sumcnd=0.0
      sumfix=0.0
      pratenew=0.0
      lpratenew=0.0
      fpratenew=0.0
      sumepr=0.0
      sumgpr=0.0
      sumqv0=0.0
      sumql0=0.0
      sumqi0=0.0
      sumqv1=0.0
      sumql1=0.0
      sumqi1=0.0

!     if(iflag.eq.1) then
!       sumqtot = 0.0
!       do k = 1,nl
!       sumqtot = sumqtot + (q(k)+clw(k)+ciw(k)+dlf(k)*delt)*dp(k)/g
!       enddo
!       write(6,*)'before loop 500, sumqtot=',sumqtot

!       sumqtot = 0.0
!       do k = 1,nl
!         sumqtot = sumqtot + (gp(k)-epl(k)-epi(k))*delt*dp(k)/g
!       end do
!       sumqtot = sumqtot - prate*delt
!       write(6,*)'sum of prec processes =',sumqtot

!       sumqtot = 0.0
!       do k = 1,nl
!         sumqtot = sumqtot + dl_ent(k)*delt*dp(k)/g
!       enddo
!       write(6,*)'sum of dl_ent =',sumqtot

!       sumqtot = 0.0
!       do k = 1,nl
!         sumqtot = sumqtot + dq_ent(k)*delt*dp(k)/g
!       enddo
!       write(6,*)'sum of dq_ent =',sumqtot
!     endif

      do 500 i=nl,1,-1

!	first update cloud amount using analytic integration
!		- sum terms for eq(37) tiedtke (1993)
!		- calculate new cloud fraction a(i)
!		- keep a(i) between 0 and 1

      aa=aa1(i)+aa2(i)+aa3(i)
      bb=bb1(i)
      if(aa+bb.eq.0.0) then
        a(i)=a(i) + aa*delt
      else if((aa+bb)*delt.gt.50.0) then
        a(i)=aa/(aa+bb)
      else
        a(i)=a(i)*exp(-(aa+bb)*delt) +		!update a(i) using analytic integration
     1                  (1.0-exp(-(aa+bb)*delt))*aa/(aa+bb) 	
      endif
      a(i)=min(1.0,a(i))
      a(i)=max(0.0,a(i))
      if(a(i).lt.1.0e-05) a(i)=0.0

!	next, update clw and ciw using finite differences
!	no precipitation yet at this point

!	partition source terms using prescribed fliq function
!	partition sink terms using existing liquid and ice ratios

      fl = fliq(t(i))
      scv_l = fl*scv(i)
      scv_i = (1.0_r8 - fl)*scv(i)
      sbl_l = fl*sbl(i)
      sbl_i = (1.0_r8 - fl)*sbl(i)
      cond_l = fl*cond(i)
      cond_i = (1.0_r8 - fl)*cond(i)
      e_l = fliqex(i)*e(i)
      e_i = ficeex(i)*e(i)
      if(dl_ent(i) .lt. 0.0) then			! here dl_ent is a sink
        dlent_l = fliqex(i)*dl_ent(i)
        dlent_i = ficeex(i)*dl_ent(i)
      else						! here dl_ent is a source
        dlent_l = fl*dl_ent(i)
        dlent_i = (1.0_r8 - fl)*dl_ent(i)
      endif

      clw(i) = clw(i) + delt*(scv_l + sbl_l + cond_l - e_l + dlent_l)
      if(clw(i).lt.0.0) then		!assumes this happens when e_l exceeds supply
        e_l = e_l + clw(i)/delt
        clw(i) = 0.0
      endif
      ciw(i) = ciw(i) + delt*(scv_i + sbl_i + cond_i - e_i + dlent_i)

      if(ciw(i).lt.0.0) then		!assumes this happens when e_i exceeds supply
        e_i = e_i + ciw(i)/delt
        ciw(i) = 0.0
      endif

      lw(i) = clw(i) + ciw(i)
      if(lw(i) .eq. 0.0) a(i) = 0.0

!	update temperature and humidity changes
!	at this point, nmo precipitatino effects yet

      dqliw(i) = delt*(-sbl_l - cond_l + e_l ) +
     1           delt*(-sbl_i - cond_i + e_i ) +
     2           delt*dq_ent(i)

      dtliw(i) = delt*((sbl_l + cond_l - e_l )*l/cpa) +
     1           delt*((sbl_i + cond_i - e_i )*(l+lf)/cpa) +
     2           delt*dt_ent(i) + 
     3           delt*scv_i*lf/cpa				! this last term is for detrained liquid water freezing

 500  continue

!	call precip2 subroutine using updated values of clw and ciw

      call precip2(pk,t,q,clw,ciw,dz,dp,rho,a,delt,dqidt_p,
     1             dqldt_p,dqdt_p,dtdt_p,lprate,fprate,paut,
     2             pcll,pacc,peva,pflx,psub,lfrc)

      prate = lprate + fprate

!	update clw and ciw due to precipitation effects
!	update change to T and q due to precipitation

      do 550 i=nl,1,-1

      clw(i) = clw(i) + delt*dqldt_p(i)
      e_l = 0.0
      if(clw(i) .lt. 0.0) then
        e_l = clw(i)/delt
        clw(i) = 0.0
      endif

      ciw(i) = ciw(i) + delt*dqidt_p(i)
      e_i = 0.0
      if(ciw(i) .lt. 0.0) then
        e_i = ciw(i)/delt
        ciw(i) = 0.0
      endif

      lw(i) = clw(i) + ciw(i)
      if(lw(i) .eq. 0.0) a(i)=0.0

      dqliw(i) = dqliw(i) + delt * ( dqdt_p(i) + e_l + e_i )
      dtliw(i) = dtliw(i) + delt * ( dtdt_p(i) - (e_l*l/cpa) - 
     1           (e_i*(l+lf)/cpa) )

 550  continue 
!*********************************************************************

      wvcol=0.0
      wlcol=0.0
      wicol=0.0
      dqcol=0.0
      dqcolabs=0.0
      do 68 i=1,nl
        dqcol = dqcol + dqliw(i)*dp(i)/g
        dqcolabs = dqcolabs + abs(dqliw(i))*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
 68   continue
      tw_end = dqcol + wvcol + wlcol + wicol + prate*delt
      dqadjr = tw_end - tw_start

!     if(dqcol .eq. 0.0) then
!       dqadjr = 0.0
!     else
!       dqadjr = (tw_end - tw_start) / dqcol
!     endif
!     if(iflag.eq.1) then
!       write(6,*)'*****************************************'
!       write(6,*)'INSIDE T93CLDS, BEFORE Q ADJUSTMENT'
!       write(6,*)'DQCOL = ', dqcol
!       write(6,*)'wvcol = ', wvcol
!       write(6,*)'wlcol = ', wlcol
!       write(6,*)'wicol = ', wicol
!       write(6,*)'precp = ', prate*delt
!       write(6,*)'tw_end= ', tw_end
!       write(6,*)'dqadjr= ', dqadjr
!     endif 
      qcol = wvcol
      wvcol=0.0
      wlcol=0.0
      wicol=0.0
      dqcol=0.0
      do 69 i=1,nl
        dqbef = dqliw(i)
        if(dqcolabs .ne. 0.0) then
          dqliw(i) = dqliw(i) - 
     1               (dqadjr*abs(dqliw(i)*dp(i)/g)/dqcolabs)*g/dp(i)
        endif
        dqcol = dqcol + dqliw(i)*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
!       if(iflag.eq.1) then
!         write(6,3820)i,dqbef,dqliw(i),q(i),q(i)+dqliw(i),clwb(i),clw(i),ciwb(i),ciw(i)
!3820     format('i=',i2,' dqbef=',1pe9.2,' dqaft=',e9.2,' qbef=',e9.2,' qaft=',e9.2,  &
!               ' clwb=',e9.2,' clwa=',e9.2,' ciwb=',e9.2,' ciwa=',e9.2)
!       endif
 69   continue
      tw_end = dqcol + wvcol + wlcol + wicol + prate*delt
!     if(iflag.eq.1) then
!       write(6,*)'INSIDE T93CLDS, AFTER Q ADJUSTMENT'
!       write(6,*)'DQCOL = ', dqcol
!       write(6,*)'wvcol = ', wvcol
!       write(6,*)'wlcol = ', wlcol
!       write(6,*)'wicol = ', wicol
!       write(6,*)'precp = ', prate*delt
!       write(6,*)'tw_end= ', tw_end
!     endif 

      secol=0.0
      dqcol=0.0
      dtcol=0.0
      dtcolabs=0.0
      do 78 i=1,nl
        dtcol = dtcol + dtliw(i)*cpa*dp(i)/g
        dtcolabs = dtcolabs + abs(dtliw(i))*cpa*dp(i)/g
        dqcol = dqcol + dqliw(i)*dp(i)/g
        secol = secol + cpa*t(i)*dp(i)/g
 78   continue
      te_end = dtcol + secol +(l+lf)*(wvcol+dqcol) + 
     1         lf*wlcol + lprate*delt*lf

      dtadjr = te_end - te_start

!     if(dtcol .eq. 0.0) then
!       dtadjr = 0.0
!     else
!       dtadjr = (te_end - te_start) / dtcol
!     endif
!     if(iflag.eq.1) then
!       write(6,*)'*****************************************'
!       write(6,*)'INSIDE T93CLDS, BEFORE T ADJUSTMENT'
!       write(6,*)'DQCOL = ', dqcol
!       write(6,*)'DTCOL = ', dtcol
!       write(6,*)'te_end= ', te_end
!       write(6,*)'dtadjr= ', dtadjr
!     endif 
      dtcol=0.0
      do 79 i=1,nl
        if(dtcolabs .ne. 0.0) then
          dtliw(i) = dtliw(i) - 
     1      (dtadjr*abs(dtliw(i)*cpa*dp(i)/g)/dtcolabs)*g/dp(i)/cpa
        endif
        dtcol = dtcol + dtliw(i)*cpa*dp(i)/g
 79   continue
      te_end = dtcol + secol +(l+lf)*(wvcol+dqcol) + 
     1         lf*wlcol + lprate*delt*lf
!     if(iflag.eq.1) then
!       write(6,*)'INSIDE T93CLDS, AFTER T ADJUSTMENT'
!       write(6,*)'DTCOL = ', dtcol
!       write(6,*)'te_end= ', te_end
!     endif 

      prec = prate*86400.0
      precl = lprate*86400.0
      precf = fprate*86400.0

      return
      end

!****************************************************************************************
      subroutine getchars(thsf  ,ttop  ,pbase ,pmidc ,ptop  ,
     1                    rlu   ,rld   ,l     ,cpa   ,gam   ,
     2                    pref  )

! 	Purpose: determine the characteristics of updrafts and downdrafts 
!		 within boundary layer clouds. This subroutine is called 
!		 when boundary layer clouds is contained within a single
!		 model layer above the boundary layer. 

      use shr_kind_mod,  only: r8 => shr_kind_r8

      implicit none

!	Input 

      real(r8), intent(in) :: thsf		! potential temperature near surface
      real(r8), intent(in) :: ttop		! temperature near top of b.l. cloud
      real(r8), intent(in) :: pbase		! pressure at base of boundary layer
      real(r8), intent(in) :: pmidc		! pressure in middle of boundary layer
      real(r8), intent(in) :: ptop		! pressure at top of boundary layer
      real(r8), intent(in) :: l			! latent heat of vaporization
      real(r8), intent(in) :: cpa		! specific heat of dry air
      real(r8), intent(in) :: gam		! R/Cp
      real(r8), intent(in) :: pref		! reference pressure

!	Output

      real(r8), intent(out) :: rlu		! liquid water content of updraft
      real(r8), intent(out) :: rld		! liquid water content of downdraft

!	Local

      real(r8) :: thbase
      real(r8) :: tbase
      real(r8) :: thcbase
      real(r8) :: tcbase
      real(r8) :: qsbase 
      real(r8) :: phalf 
      real(r8) :: thcg
      real(r8) :: thovertmid
      real(r8) :: tcg
      real(r8) :: qsg
      real(r8) :: tcmid
      real(r8) :: qsmid
      real(r8) :: dqs
      real(r8) :: dtc
      real(r8) :: thcnew
      real(r8) :: thcld
      real(r8) :: dt
      real(r8) :: thtop
      real(r8) :: thctop 
      real(r8) :: tctop
      real(r8) :: qstop
      real(r8) :: qsat
      integer :: it

!	first calculate the updraft characteristics

!	surface air rises dry adiab. to base of b.l. cloud

      thbase=thsf			!pot. temperature of rising parcel at base of b.l. cloud
      tbase=thbase*(pbase/pref)**gam	!temperature of rising parcel at base of b.l. cloud

!	once inside cloud, parcel now rises on moist adiabat

      thcbase=thbase			!cloud pot. temperature at base of b.l. cloud
      tcbase=tbase			!cloud temperature at base of b.l. cloud
      qsbase=qsat(tcbase,pbase)		!saturation humidity of parcel at base of bl cloud

!	itertate to find potential temp on moist adiabat at pres=pmidc

      it=0				!initialize iteration counter
      phalf=0.5*(pmidc+pbase)		!pressure halfway between cloud base and mid cloud height
      thcg=thcbase+1.0			!initial guess for cloud potential temp.
      thovertmid=((pref/phalf)**gam)
 111  it=it+1
      tcg=thcg*((pmidc/pref)**gam)
      qsg=qsat(tcg,pmidc)
      tcmid=0.5*(tcbase+tcg)
      qsmid=0.5*(qsbase+qsg)
      dqs=qsg-qsbase
      dtc=tcg-tcbase
      thcnew=thcbase-thovertmid*l*(dqs-(qsmid*dtc/tcmid))/cpa
      if(abs(thcnew-thcg).gt.0.01.and.it.lt.100) then
        thcg=0.25*thcnew+0.75*thcg
        goto 111
      endif
      thcld=0.5*(thcnew+thcg)		!cloud potential temperature at p=pmidc
      dt=thcld-thbase			!temperature diff between base and midpoint
      rlu=dt*cpa/l			!condensed water necessary for temperature rise dt

!	now calculate the downdraft characteristics

      thtop=ttop*((pref/ptop)**gam)	!pot. temperature of sinking parcel at top of b.l. cloud

!	air just above bl cloud sinks along moist adiabat to pres=pmidc

      thctop=thtop			!cloud pot. temperature at top of b.l. cloud
      tctop=ttop			!cloud temperature at top of b.l. cloud
      qstop=qsat(tctop,ptop)		!saturation humidity of parcel at top of bl cloud

!	itertate to find potential temp on moist adiabat at pres=pmidc

      it=0				!initialize iteration counter
      phalf=0.5*(pmidc+ptop)		!pressure halfway between cloud top and mid cloud height
      thcg=thctop-1.0			!initial guess for cloud potential temp.
      thovertmid=((pref/phalf)**gam)
 222  it=it+1
      tcg=thcg*((pmidc/pref)**gam)
      qsg=qsat(tcg,pmidc)
      tcmid=0.5*(tctop+tcg)
      qsmid=0.5*(qstop+qsg)
      dqs=qsg-qstop
      dtc=tcg-tctop
      thcnew=thctop-thovertmid*l*(dqs-(qsmid*dtc/tcmid))/cpa
      if(abs(thcnew-thcg).gt.0.01.and.it.lt.100) then
        thcg=0.25*thcnew+0.75*thcg
        goto 222
      endif
      thcld=0.5*(thcnew+thcg)		!cloud potential temperature at p=pmidc
      dt=thcld-thtop			!temperature diff between top and midpoint
      rld=dt*cpa/l			!condensed water necessary for temperature rise dt

      return
      end

!-------------------------------------------------------------
      subroutine getthc2s(p1   ,p2   ,th1  ,th2  ,l    ,
     1                    cpa  ,gam  ,pref )

!	itertate to find potential temp at point 2 on 
!	moist adiabat connecting pressures p1 and p2

      use shr_kind_mod,  only: r8=>shr_kind_r8

      implicit none

!	Input

      real(r8), intent(in) :: p1
      real(r8), intent(in) :: p2
      real(r8), intent(in) :: th1
      real(r8), intent(in) :: l
      real(r8), intent(in) :: cpa
      real(r8), intent(in) :: gam
      real(r8), intent(in) :: pref

!	Output

      real(r8), intent(out) :: th2

!	Local

      real(r8) :: phalf
      real(r8) :: thcg
      real(r8) :: thovertmid
      real(r8) :: t1
      real(r8) :: qs1
      real(r8) :: tcg
      real(r8) :: qsg
      real(r8) :: tcmid
      real(r8) :: qsmid
      real(r8) :: dqs
      real(r8) :: dtc
      real(r8) :: thcnew
      real(r8) :: qsat			! real(r8) :: function

      integer :: it

      it=0				!initialize iteration counter
      phalf=0.5*(p2+p1)			!pressure halfway between cloud base and mid cloud height
      thcg=th1+1.0			!initial guess for cloud potential temp.
      thovertmid=((pref/phalf)**gam)
      t1=th1*(p1/pref)**gam
      qs1=qsat(t1,p1)

 111  it=it+1
      tcg=thcg*((p2/pref)**gam)
      qsg=qsat(tcg,p2)
      tcmid=0.5*(th1+tcg)
      qsmid=0.5*(qs1+qsg)
      dqs=qsg-qs1
      dtc=tcg-t1
      thcnew=th1-thovertmid*l*(dqs-(qsmid*dtc/tcmid))/cpa
      if(abs(thcnew-thcg).gt.0.01.and.it.lt.100) then
        thcg=0.25*thcnew+0.75*thcg
        goto 111
      endif
      th2=0.5*(thcnew+thcg)		!cloud potential temperature at p=pmidc

      return
      end
!************************************************************************************
      subroutine precip2(pk,t,q,ql,qi,dz,dp,rho,cld,delt,
     1                  dqidt,dqldt,dqdt,dtdt,           
     2                  rlsurf,rfsurf,paut,pcll,pacc,peva,
     3                  pflx,psub,lfrc)

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use ppgrid,        only: nl => pver
      use physconst,     only: l => latvap,  g => gravit, cpa =>cpair,
     1                         lf => latice 
      implicit none

!	input quantities

      real(r8), intent(in) :: q(nl)			! specific humidity (kg/kg)
      real(r8), intent(in) :: ql(nl)			! cloud liquid water mixing ratio (kg/kg)
      real(r8), intent(in) :: qi(nl)			! cloud ice water mixing ratio (kg/kg)
      real(r8), intent(in) :: t(nl)			! temperature (deg K)
      real(r8), intent(in) :: pk(nl)			! mid-layer pressure (Pa)
      real(r8), intent(in) :: rho(nl)			! air density (kg/m**3)
      real(r8), intent(in) :: dp(nl)			! layer pressure thickness (Pa)
      real(r8), intent(in) :: dz(nl)			! layer geometric thickness (m)
      real(r8), intent(in) :: cld(nl)			! cloud fraction (fraction)
      real(r8), intent(in) :: delt
      real(r8), intent(in) :: lfrc			! land fraction

!	output quantities

      real(r8), intent(out) :: dqidt(nl)		! change in cloud ice water due to precip processes (kg/kg/sec)
      real(r8), intent(out) :: dqldt(nl)		! change in cloud liquid water due to precip processes (kg/kg/sec)
      real(r8), intent(out) :: dtdt(nl)			! change in temp due to precipitation phase changes
      real(r8), intent(out) :: dqdt(nl)			! change in humidity due to precip processes (kg/kg/sec)
      real(r8), intent(out) :: rlsurf			! liquid precipitation rate at surface (kg/m**2/sec)
      real(r8), intent(out) :: rfsurf			! frozen precipitation rate at surface (kg/m**2/sec)
      real(r8), intent(out) :: paut(nl)		! conversion of cloud water to precip due to autoconversion (g/m**2/sec)
      real(r8), intent(out) :: pcll(nl)		! conversion of cloud water to precip due to collection (g/m**2/sec)
      real(r8), intent(out) :: pacc(nl)		! conversion of cloud water to precip due to accretion (g/m**2/sec)
      real(r8), intent(out) :: peva(nl)		! evaporation of liquid precipitation (g/m**2/sec)
      real(r8), intent(out) :: pflx(nl)		! flux divergence of ice particles    (g/m**2/sec)
      real(r8), intent(out) :: psub(nl)		! sublimation of frozen precipitation (g/m**2/sec)

!	local quantities
 
      real(r8) :: ls					! sum of latent heat of sublimation (l + lf)
      real(r8) :: rfmelt				! amount of frozen precip that is melting in a layer (kg/m**2/sec)
      real(r8) :: rfmeltmax				! maximum value of rfmelt allowed in a layer
      real(r8) :: parea					! fractional area in which precipitation is falling
      real(r8) :: parea1
      real(r8) :: vf					! fall speed of frozen precipitation (m/sec)
      real(r8) :: qs					! saturation specific humidity for a given layer
      real(r8) :: es					! saturation vapor pressure for a given layer
      real(r8) :: qs1					! used to calculate dqsdt
      real(r8) :: qs2					! used to calculate dqsdt
      real(r8) :: dqsdt					! derivative of qs w.r.t. temperature for a given layer
      real(r8) :: rf					! grid-mean snowfall rate (kg/m**2/sec)
      real(r8) :: rfloc					! local snowfall rate (kg/m**2/sec)
      real(r8) :: app					! term used in calculating sublimation of falling ice
      real(r8) :: bpp					! term used in calculating sublimation of falling ice
      real(r8) :: dev					! term used in calculating sublimation of falling ice
      real(r8) :: fev					! term used in calculating sublimation of falling ice
      real(r8) :: qnew					! term used to temporarily hold new value of q
      real(r8) :: epi					! sublimation of falling ice (snow) (kg/kg)
      real(r8) :: epimax				! maximum allowed sublimation in a layer
      real(r8) :: rl					! grid-mean rainfall rate (kg/m**2/sec)
      real(r8) :: rlloc					! local rainfall rate (kg/m**2/sec)
      real(r8) :: acorain				! rainfall rate in layer due to autoconversion (kg/m**2/sec)
      real(r8) :: aprime				! term used in calculating evaporation of falling rain
      real(r8) :: bprime				! term used in calculating evaporation of falling rain
      real(r8) :: bev					! term used in calculating evaporation of falling rain
      real(r8) :: cev					! term used in calculating evaporation of falling rain
      real(r8) :: epl					! evaporation of falling rain (kg/kg)
      real(r8) :: eplmax				! maximum allowed evaporation in a layer
      real(r8) :: alpha					! term used in calculating precipitation of cloud ice
      real(r8) :: r1   					! term used in calculating precipitation of cloud ice
      real(r8) :: r2   					! term used in calculating precipitation of cloud ice
      real(r8) :: rffall				! precipitation entering top of layer that falls through bottom of layer (kg/kg)
      real(r8) :: rfstay				! precipitation entering top of layer that remains in layer (kg/kg)
      real(r8) :: qifall				! ice initially in a layer that falls out (kg/kg)
      real(r8) :: aco					! auto-conversion of cloud liquid water to precip (kg/kg) 
      real(r8) :: qlc					! in-cloud value of cloud liquid water
      real(r8) :: cau					! term used in calculating auto-conversion of cloud liquid water
      real(r8) :: qcr					! critical mixing ratio of q used in auto-conversion calculation
      real(r8) :: qval					! temporary value of q (auto-conversion)
      real(r8) :: sfac					! slope factor in size distribution of falling ice (accretion)
      real(r8) :: accr					! accretion rate of liquid water to ice (1/sec)
      real(r8) :: acc					! change in ql due to accretion (kg/kg)
      real(r8) :: colr					! rate of collection of cloud liquid water by liquid precip (1/sec)
      real(r8) :: colr1					! rate of collection of cloud liquid water by liquid precip (1/sec)
      real(r8) :: colr2					! rate of collection of cloud liquid water by liquid precip (1/sec)
      real(r8) :: col					! change in ql due to collection (kg/kg)
      real(r8) :: qlnew					! temporary value of ql
      real(r8) :: ratio					! used if necessary to reduce cloud liquid water sinks
      real(r8) :: rfbot					! falling ice precipitation through bottom of layer (kg/kg)
      real(r8) :: rlbot					! falling liquid precipitation through bottom of layer (kg/kg)
      real(r8) :: qsat					! function
      real(r8) :: esat					! function
      real(r8) :: dropcn				! droplet concentration (m**-3) (derived using lfrc)
      integer :: k

!	constants and parameters

      real, parameter :: pi = 3.14159
      real, parameter :: ka = 2.40e-02			! thermal conductivity (J/m/s/K)
      real, parameter :: rv = 461.0			! specific gas constant for water vapor (J/kg/K)
      real, parameter :: chi = 2.21e-05			! diffusivity of water vapor in air at 1000hPa (m**2/sec)
      real, parameter :: eac = 0.70			! falling ice / cloud liquid water collection efficiency (accretion)
      real, parameter :: eco = 0.70			! falling rain / cloud liquid water collection efficiency (collection)
      real, parameter :: rhof = 100.0			! bulk density of falling ice (kg/m**3)
      real, parameter :: rmin = 1.0e-12			! minimum precipitation value (kg/m**2/sec)
      real, parameter :: dropcn_l = 2.0e+08		! cloud droplet conc over land (m**-3)
      real, parameter :: dropcn_o = 7.5e+07		! cloud droplet conc over water (m**-3)
      real, parameter :: rhowat = 1000.0		! density of water (kg/m**3)
      real, parameter :: visc = 1.72e-05		! viscosity of air (kg/m/sec)
      real, parameter :: ec = 0.55			! mean collector efficiency for autoconversion
      real, parameter :: rcr = 7.5e-06			! critical radius (meters)

      ls = l + lf
      dropcn = dropcn_l*lfrc + dropcn_o*(1.0_r8 - lfrc)

      vf=0.0							!initialize fall speed of frozen precip from above (m/sec)
      rf=0.0							!initialize frozen precipitation rate from above layer (kg/m**2/sec)
      rl=0.0							!initialize liquid precipitation rate from above layer (kg/m**2/sec)
      parea=0.0							!initialize fractional area in which precipitation is falling
      do 3 k=nl,1,-1

      dtdt(k) = 0.0						!initialize temperature change due to phase change of precipitation
      qs1 = qsat(t(k)+0.25,pk(k))
      qs2 = qsat(t(k)-0.25,pk(k))
      dqsdt=(qs1-qs2)/0.5
      qs = qsat(t(k),pk(k))
      es = esat(t(k),pk(k))

!	check if any falling ice entering layer with t>275.15
!	if so, melt ice precip, but do not allow layer temp to
!	increase above 275.15

      if(rf .gt. 0.0 .and. t(k) .gt. 275.15) then
        rfmeltmax = (t(k)-275.15)*cpa*(dp(k)/g)/lf/delt		!maximum allowed melting (kg/m**2/sec) without letting temp go below 275.15
        rfmelt = min(rf, rfmeltmax)
        rf = rf - rfmelt
        rl = rl + rfmelt
        dtdt(k) = dtdt(k) - rfmelt*lf/cpa/(dp(k)/g)		!cooling due to melting (K/sec)
      endif

!	check to see if parea is undefined if there is any precipitating
!	ice or liquid. this may occur if there is a cloud-free layer 
!	containing some cloud ice without any precip falling into
!	the layer. if this occurs, use a value of 0.5 for parea. 
!	this is crude approximation, but should not be
!	too important as these types of layers should not 
!	contain too much cloud ice.

      parea1 = parea
      if(rl+rf .gt. 0.0 .and. parea .eq. 0.0) parea1 = 0.5

!	calculate sublimation of falling ice
!	do this only if there is ice falling into layer (rf>0)
!	use fall speed (vf) from above layer
!	update rf

      rfloc = 0.0
      if(rf .gt. 0.0) rfloc = rf/parea1			!local snowfall (ice) rate (kg/m**2/sec)
      if(rfloc .gt. 0.0) then
        app = (ls/ka/t(k))*((ls/rv/t(k)) - 1.0)
        bpp = rv*t(k)/chi/es
        dev = (4.0*rf/(rho(k)*(app+bpp)*qs*pi*rhof*vf))* 		!here no parea factor since it cancels with rf raised to power of 1
     1   ((0.65*sfac*sfac)+0.493*sfac*sqrt(sfac)*sqrt(vf*rho(k)/visc))
        fev = 1.0 + 0.5*dev*delt*(1.0 + (ls/cpa)*dqsdt)
        qnew = q(k) + (dev*delt/fev)*(qs-q(k)) 
        qnew = min(qnew, qs)					!insure new humidity not supersaturated
        epi = qnew - q(k)					!sublimation of snow (kg/kg) (should be a positive number)
        epimax = (1.0-cld(k))*rf*delt/(dp(k)/g)			!this is max epi since sublimation only allowed in cloud free area
        epi = min(epi,epimax)
      else
        epi = 0.0
      endif

      rf = rf - epi*(dp(k)/g)/delt				!update rf. this value is now used for other processes in layer (kg/m**2/sec)
      if(rf .lt. rmin) then					!check if this is above minimum value
        epi = epi + rf*delt/(dp(k)/g)				!sublimate this small precip amount (kg/kg)
        rf = 0.0
      endif

!	calculate evaporation of rain
!	do this only if there is rain falling into layer (rl>0)
!	update rl

      rlloc = 0.0
      if(rl .gt. 0.0) rlloc = rl/parea1				!local rainfall rate kg/m**2/sec
      if(rlloc .gt. 0.0) then
        aprime=(l/ka/t(k))*((l/rv/t(k)) - 1.0)
        bprime=rv*t(k)/chi/es
        cev = parea1*870.0*(rlloc**0.61)/ 
     1                   (sqrt(rho(k))*(aprime+bprime)*qs)
        bev = 1.0 + 0.5*cev*delt*(1.0 + (l/cpa)*dqsdt)
        qnew = q(k) + (cev*delt/bev)*(qs-q(k))
        qnew = min (qnew, qs)					!insure new humidity not supersaturated
        epl = qnew - q(k)					!evaporation of rain (kg/kg) (should be a positive number)
        eplmax = (1.0-cld(k))*rl*delt/(dp(k)/g)			!this is max epl since evaporation only allowed in cloud free areas
        epl = min(epl,eplmax)
      else
        epl = 0.0
      endif

      rl = rl - epl*(dp(k)/g)/delt				!update rl. this updated value is now used for other processes in layer
      if(rl .lt. rmin) then					!check if this is above minimum value
        epl = epl + rl*delt/(dp(k)/g)				!evaporate this small precip amount (kg/kg)
        rl = 0.0
      endif

!	calculate precipitation of cloud ice from layer and
!	the fraction of falling ice from above that falls 
!	through layer.
!
!       update vf if any cloud ice in layer. if no cloud ice, the
!	use vf from above layer. also use vf from above layer if
!	there is no cloud in this layer. this would occur if cloud
!	ice falls into and remains in a cloud-free layer.	

      if(qi(k).gt.0.0.and.cld(k).ne.0.0) then
        vf = 3.23*(rho(k)*qi(k)/cld(k))**0.17			!fall speed of frozen precip (m/sec)
      endif
      if(vf .gt. 0.0) then
        alpha = -vf*delt/dz(k)					!switch sign since dz is negative here
        r1 = exp(-alpha)
        r2 = (1.0 - exp(-alpha))/alpha
      else							!here no ice in layer and no ice falling 
        r1 = 1.0						!no generation of frozen precip from this layer
        r2 = 1.0						!set r2=1, but will not matter as rf should equal 0
      endif
      rffall = rf*delt*(1.0-r2)/(dp(k)/g)			!precip entering top of layer that falls through bottom of layer (kg/kg)
      rfstay = rf*delt*r2/(dp(k)/g)				!precip entering top of layer that remains in layer (kg/kg)
      qifall = qi(k)*(1.0-r1)					!ice initially present in layer that falls out (kg/kg)

!	calculate autoconversion of cloud liquid water
!	use manton-cotton formulation
!	use rotstayn equations (mks units)

      qlc=0.0
      if(cld(k).gt.0.0) qlc=ql(k)/cld(k)
      if(qlc .eq. 0.0) then
        aco = 0.0
      else
        cau = (0.104*9.81*ec*rho(k)**(4./3.))/
     1            (visc*(dropcn*rhowat)**(1./3.))
        qcr = 4.0*3.14159*rhowat*dropcn*(rcr**3.0)/(3.0*rho(k))	!critical mixing ratio (kg/kg)
        if(qlc .le. qcr) then
          aco = 0.0
        else
          qval = ( (qlc**(-4.0/3.0)) + 4.0*cau*delt/3.0 ) **(-0.75)
          qnew = cld(k)*max(qcr, qval)				!grid-mean value; units=kg/kg 
          aco = qnew - ql(k)					!analytic evaluation of autoconversion; should be negative value 
          aco = -aco						!make a positive value (units=kg/kg)
        endif
      endif

!	calculate accretion
!	make sure accretion does not exceed 50% of available cloud liquid water

      sfac = 1.6e+03 * 10.0**(0.023*(273.15-t(k)))		!slope factor in size dist. of falling ice
      accr = eac*sfac*rf/(2.0*rhof)				!rate of accretion of liq water to ice (1/sec)
!     if(t(k) .gt. 273.15) accr = 0.0				!set to 0 if temp above melt point
      qnew = ql(k)*(1.0-0.5*accr*delt)/(1.0+0.5*accr*delt)	!new cloud liq water due only to accretion; center-weighted difference
      acc = qnew - ql(k)					!change in ql due to accretion (should be negative number)
      acc = -acc						!change to positive number (kg/kg)
      if(acc .gt. 0.5*ql(k)) acc=0.5*ql(k)			!limit accretion to 50% of total cloud liquid water 

!	calculate collection 

      rlloc = 0.0
      if(rl .gt. 0.0) rlloc = rl/parea1				!local rainfall rate kg/m**2/sec
      colr1 = 0.24*parea1*eco*(rlloc**0.75)			!rate of collection of cloud liquid water by liquid precip (1/sec)
      acorain = aco*dp(k)/g/delt				!rate of rainfall produced in this layer due to autoconversion (kg/m**2/sec)
      colr2 = 0.5*0.24*1.00*eco*(acorain**0.75)			!collection due to rain produced in this layer (0.5 for midlayer)
      colr = colr1 + colr2
      qnew = ql(k)*(1.0-0.5*colr*delt)/(1.0+0.5*colr*delt)	!new cloud liq water due to only collection; center-weighted difference
      col = qnew - ql(k)					!change in ql due to collection (should be negative number)
      col = -col						!change to positive number (kg/kg)
      if(col .gt. 0.5*ql(k)) col=0.5*ql(k)			!limit collection to 50% of total cloud liquid water

!	sum terms for this layer

      dqidt(k) = -(qifall - rfstay) / delt			!change in cloud ice water due to precip processes (kg/kg/sec)
      pflx(k) = -1000.0*dqidt(k)*dp(k)/g			!flux divergence of ice in this layer (g/m**2/sec)

 341  dqldt(k) = -(acc+col+aco) / delt				!change in cloud liquid water due to precip processes (kg/kg/sec)
      qlnew = ql(k) + delt*dqldt(k)
      if(qlnew .lt. 0.0) then					!here need to reduce cloud liquid water sinks
        ratio=ql(k)/(ql(k)-qlnew)				!calculate reduction ratio
        acc = ratio*acc						!and apply to all terms evenly
        col = ratio*col
        aco = ratio*aco
        qlnew = 0.0
      endif
      paut(k) = 1000.0*aco*dp(k)/delt/g				!compute conversion term for diagnostics (g/m**2/sec)
      pcll(k) = 1000.0*col*dp(k)/delt/g				!compute conversion term for diagnostics (g/m**2/sec)
      pacc(k) = 1000.0*acc*dp(k)/delt/g				!compute conversion term for diagnostics (g/m**2/sec)

      rfbot = qifall + rffall + acc 			 	!falling frozen precip through bottom of layer; units = kg/kg
      rfbot = rfbot*(dp(k)/g)/delt				!change units = kg/m**2/sec
      if(rfbot .lt. rmin) then					!check if this is above minimum value
        epi = epi + rfbot*delt/(dp(k)/g)			!sublimate this small precip amount (kg/kg)
        rfbot = 0.0
      endif

      rlbot = rl + (col+aco)*(dp(k)/g)/delt			!falling liquid precip through bottom of layer; units = kg/m**2/sec
      if(rlbot .lt. rmin) then					!check if this is above minimum value
        epl = epl + rlbot*delt/(dp(k)/g)			!evaporate this small precip amount (kg/kg)
        rlbot = 0.0
      endif
      peva(k) = 1000.0*epl*dp(k)/delt/g				!compute conversion term for diagnostics (g/m**2/sec)
      psub(k) = 1000.0*epi*dp(k)/delt/g				!compute conversion term for diagnostics (g/m**2/sec)

      dtdt(k) = dtdt(k) +        
     1          acc*lf/cpa/delt - 				!dT (warming) due to accretion (liquid cloud water to frozen precip) (K/sec)
     2          epl*l/cpa/delt -  				!dT (cooling) due to evaporation (liquid cloud water to vapor) (K/sec) 
     3          epi*ls/cpa/delt					!dT (cooling) due to sublimation (frozen cloud water to vapor) (K/sec) 

      rf = rfbot
      rl = rlbot
      if(rf+rl .gt. 0.0) then
        parea = max(parea, cld(k))				!precipitating fraction of gridbox for precip falling from this layer 
      else
        parea = 0.0
      endif

      dqdt(k) = (epi/delt) + (epl/delt)				!rate of change of q due to sublimation and evaporation of precip (kg/kg/sec)

 3    continue

      rlsurf = rl						!liquid precipitation rate at surface (kg/m**2/sec)
      rfsurf = rf						!frozen precipitation rate at surface (kg/m**2/sec)
 
      return
      end 

!*******************************************************************
      subroutine lcl(pk,th,q,lev,pval,pref)

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use ppgrid,        only: nl => pver

      implicit none

!	determine the lifting condensation level using iterative technique.  

!	Input

      real(r8), intent(in) :: pk(nl)
      real(r8), intent(in) :: th(nl)
      real(r8), intent(in) :: q(nl)
      real(r8), intent(in) :: pref
      integer, intent(in) :: lev

!	Output

      real(r8), intent(out) :: pval

!	Local

      real(r8) :: pg		! guess value of plcl
      real(r8) :: pgold		! previous guess value
      real(r8) :: pgmb		! guess value in mb
      real(r8) :: es
      real(r8) :: lnterm
      real(r8) :: tg
      real(r8) :: diffpg
      integer :: n

      pg=pk(lev)-20000.0
 300  do 400 n=1,500
      pgmb=pg/100.0
      pgold=pg
      es=(pgmb*q(lev))/(0.622+(0.378*q(lev)))
      lnterm=log(es/6.11)
      tg=(4714.742-(35.86*lnterm))/(17.26-lnterm)
      pg=((tg/th(lev))**3.497)*pref
      diffpg=pg-pgold
      if(abs(diffpg).lt.50.0) goto 1111
 400  continue
 1111 pval=pg
      if(pval.gt.pk(lev)) pval=pk(lev)
      return 
      end


!****************************************************************************
      function esat(t,pk) result (esat_result)

!	this function returns a value for the saturated specific
!	humidity using the routine from ccm3.

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use wv_saturation, only: aqsat
      implicit none

      real(r8) :: esat_result
      real(r8), intent(in) :: t
      real(r8), intent(in) :: pk
      real(r8) :: es
      real(r8) :: qs

      integer, parameter :: ii = 1
      integer, parameter :: kk = 1
      integer, parameter :: ilen = 1
      integer, parameter :: kstart = 1
      integer, parameter :: kend = 1

      real(r8) :: tarr(ii,kk)
      real(r8) :: pkarr(ii,kk)

      real(r8) :: esarr(ii,kk)         ! Saturation vapor pressure
      real(r8) :: qsarr(ii,kk)         ! Saturation specific humidity

      tarr(1,1) = t
      pkarr(1,1) = pk
      call aqsat(tarr  ,pkarr ,esarr  ,qsarr  ,ii     ,
     1           kk    ,ilen  ,kstart ,kend   )
      esat_result = esarr(ii,kk)

      return
      end


!****************************************************************************
      function qsat(t,pk) result (qsat_result)

!	this function returns a value for the saturated specific
!	humidity using the routine from ccm3.

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use wv_saturation, only: aqsat
      implicit none

      real(r8) :: qsat_result
      real(r8), intent(in) :: t
      real(r8), intent(in) :: pk
      real(r8) :: es
      real(r8) :: qs

      integer, parameter :: ii = 1
      integer, parameter :: kk = 1
      integer, parameter :: ilen = 1
      integer, parameter :: kstart = 1
      integer, parameter :: kend = 1

      real(r8) :: tarr(ii,kk)
      real(r8) :: pkarr(ii,kk)

      real(r8) :: esarr(ii,kk)         ! Saturation vapor pressure
      real(r8) :: qsarr(ii,kk)         ! Saturation specific humidity

      tarr(1,1) = t
      pkarr(1,1) = pk
      call aqsat(tarr  ,pkarr ,esarr  ,qsarr  ,ii     ,
     1           kk    ,ilen  ,kstart ,kend   )
      qsat_result = qsarr(ii,kk)

      return
      end


!****************************************************************************
      function fliq(tk) result (fliq_result)

!       this function calculates the fraction of liquid cloud water
!       1-fliq = fice, the fraction of cloud ice. taken from
!       smith, 1990

      use shr_kind_mod,  only: r8 => shr_kind_r8

      implicit none

      real(r8) :: fliq_result
      real(r8), intent(in) :: tk
      real(r8) :: tc

      tc=tk-273.15
      if(tc.ge.0.0) fliq_result=1.0
      if(tc.lt.0.0.and.tc.ge.-5.0) fliq_result=1.0-((0.2*tc)**2)/3.0
      if(tc.lt.-5.0.and.tc.gt.-15.0) 
     1   fliq_result=(((tc+15.0)/5.0)**2.0)/6.0
      if(tc.le.-15.0) fliq_result=0.0

      return
      end

