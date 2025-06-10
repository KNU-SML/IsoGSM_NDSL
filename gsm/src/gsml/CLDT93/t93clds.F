!----------------------------------------------------------------------
! original file is t93clds.F90 from ECPC SCM by S.Iacobellis
!
!   may. 2006     modified by A.Shimpo
!     * "remodified from original file with using module"
!       -> closer to original file
!     * convert freeform format to fixed format
!     * [bugfix] dqfix=0.0 initially
!     * setting "define TUNE1" : same as aug.2005 modification
!       ** add calculation of tropopause level = itrop
!       ** modify rhcrit when sgk(i)<=0.65.
!          rhcrit=0.80 -> 0.80+0.20*(0.65-sgk(i))/(0.65-sgk(itrop))
!       ** dqsdt=dqsdpma*(wpa(i)+0.05*g*cmf(i)) + ...
!         ->dqsdt=dqsdpma*(wpa(i)+g*cmf(i)) + ...
!       ** if((sac(i).gt.sacmx).or.(q(i).ge.qs)) then ...
!         ->if((sac(i).gt.sacmx).or.((dqsdt.lt.0.0).and.(q(i).ge.qs))) then ...
!       ** modify e(i) when lw(i)-e(i)*delt < 0
!       ** modify dcoef
!
!   mar. 2008     modified by A.Shimpo
!     * setting "define TUNE2". default setting for ECPC G-RSM
!
!   
!----------------------------------------------------------------------
!**********************************************************************
      subroutine t93clds(sgk   ,p     ,pk    ,t     ,th    ,
     1                   q     ,lw    ,clw   ,ciw   ,rho   ,dlf   ,du ,
     2                   z     ,dtsur ,dqsur ,solin ,wpa   ,
     3                   cmf   ,hrtr  ,hrsl  ,dthor ,hpbl  ,
     4                   fqs   ,fss   ,pref  ,delt  ,a     ,
     5                   dtliw ,dqliw ,prec  ,precl ,precf ,iflag,
     6                   sacv  ,sac   ,sabl  ,disp )

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use ppgrid,        only: nl => pver, nl1 => pverp
      use physconst,     only: l => latvap,  g => gravit, cpa => cpair,
     1                         gam => cappa, ra => rair,  lf => latice
      implicit none
 
!!#define TUNE1  ! tuning in aug, 2005

#define TUNE2 ! default for ECPC G-RSM

!	cloud and liquid water parameterization of tiedtke, 1993
!
!	written to be used in CAM3 and is called once for each column.
!	in this subroutine, layer number increases from the ground up, so
!	cam3 arrays must be inverted prior to call. 

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

!	input/output variables

      real(r8), intent(inout) :: a(nl)		! cloud fraction
      real(r8), intent(inout) :: lw(nl)		! cloud water/ice kg/kg (grid average)
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
      real(r8) scv(nl)		! cloud water source due to convection (kg/kg/sec)
      real(r8) scv_l, scv_i
      real(r8) sbl(nl)		! cloud water source due to stratocumulus (boundary layer) clouds
      real(r8) sbl_l, sbl_i
      real(r8) sacv(nl)		! cloud source due to convection (fraction/sec)
      real(r8) sac(nl)		! cloud source due to stratiform cloud (fraction/sec)
      real(r8) sabl(nl)		! cloud source due to stratocumulus (boundary layer) clouds (fraction/sec)
      real(r8) cond(nl)		! condensation rate in stratiform clouds
      real(r8) cond_l, cond_i
      real(r8) gp(nl)		! conversion of cloud (total)water to precipitation (kg/kg/sec)
      real(r8) gpl(nl)		! conversion of cloud liquid water to precipitation (kg/kg/sec)
      real(r8) gpi(nl)		! conversion of cloud ice water to precipitation (kg/kg/sec)
      real(r8) ep(nl)		! evaporation of precipitation falling into layer (kg/kg/sec)
      real(r8) epl(nl)		! evaporation of liquid precipitation falling into layer (kg/kg/sec)
      real(r8) epi(nl)		! evaporation of frozen precipitation falling into layer (kg/kg/sec)
      real(r8) e(nl)		! evaporation rate of cloud water
      real(r8) e_l, e_i
      real(r8) disp(nl)		! cloud disipation rate (fraction/sec)
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

#ifdef TUNE1
      integer :: itrop
      real(r8) :: dthdp
#endif

#ifndef TUNE1
      real, parameter :: dcoef = 1.0e-06	! diffusion coefficient for evaporation of cloud water from turbulent mixing (1/sec)
#else
      real, parameter :: dcoef = 2.0e-06	! diffusion coefficient for evaporation of cloud water from turbulent mixing (1/sec)
#endif

#ifdef TUNE2
      real(r8) :: tt1
      real(r8) :: tt2
      real(r8) dpg(nl)		! layer thickness (Pa) / g
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
      dqfix=0.0

c###test
cc      return
c###test

      do 5 i=1,nl
      if (clw(i)+ciw(i) .eq. 0.0) then
        fliqex(i) = 0.0
        ficeex(i) = 0.0
      else
        fliqex(i) = clw(i) / (clw(i) + ciw(i))
        ficeex(i) = ciw(i) / (clw(i) + ciw(i))
      endif
      dp(i)=p(i)-p(i+1)
      dz(i)=z(i)-z(i+1)			!caution: negative values here
#ifdef TUNE2
      dpg(i)=dp(i)/g
#endif
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
#ifdef TUNE2
        tt1=dpg(i)
        secol = secol + cpa*t(i)*tt1
        wvcol = wvcol + q(i)*tt1
        wlcol = wlcol + clw(i)*tt1
        wicol = wicol + ciw(i)*tt1
        dlfcol = dlfcol + dlf(i)*delt*tt1
#else
        secol = secol + cpa*t(i)*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
        dlfcol = dlfcol + dlf(i)*delt*dp(i)/g
#endif
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
#ifdef TUNE2
        exit
#else
        goto 21
#endif
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
#ifdef TUNE2
        exit                                            !either the CBL is non-existant or the cloud
#else
        goto 3839					!either the CBL is non-existant or the cloud
#endif
      endif						!layer extends across more than one model layer

      thes=th(k)*exp(l*qsat(t(k),pk(k))/(cpa*t(k)))	!saturation equivalent potential temperature of layer k

#ifdef TUNE2
      if(istart.eq.0) then
        if(thes.le.the2) then		!enter this block if we are encountering positive
          istart=1			!buoyancy (thes<the2) for the first time.
          layb=k			!set layb as the bottom of this cloud
        else
          nlay=k-layb			!the cloud layer defined as (thes>the2). set nlay 
          exit
        endif
      endif
#else
      if(istart.eq.0.and.thes.le.the2) then		!enter this block if we are encountering positive
        istart=1					!buoyancy (thes<the2) for the first time.
        layb=k						!set layb as the bottom of this cloud
      endif

      if(istart.eq.1.and.thes.gt.the2) then		!enter this block if we have reached the top of
        nlay=k-layb					!the cloud layer defined as (thes>the2). set nlay 
        goto 3839					!to the number of model layers within this cloud
      endif						!and exit loop 3838
#endif

 3838 continue
 3839 if(nlay.gt.1) goto 401				!either no condensation or cloud spans more than one model layer
      lbase=layb
      if(sacv(lbase).gt.0.0.or.scv(lbase).gt.0.0) goto 401	!skip section if convection present 
#ifdef TUNE2
      fqbase=dqsur(lbase)*dpg(lbase)				!approximate moisture flux thru base of cloud
#else
      fqbase=dqsur(lbase)*dp(lbase)/g				!approximate moisture flux thru base of cloud
#endif

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
#ifdef TUNE2
      tflxnet=dtsur(lay)*cpa*dpg(lay)		!net flux needed for layer
      qflxnet=dqsur(lay)*l*dpg(lay)		!net flux needed for layer
#else
      tflxnet=dtsur(lay)*cpa*dp(lay)/g		!net flux needed for layer
      qflxnet=dqsur(lay)*l*dp(lay)/g		!net flux needed for layer
#endif
      fs(lay+1)=fs(lay)-tflxnet
      fq(lay+1)=fq(lay)-qflxnet
 30   continue

      do 31 lev=1,lbase+1
#ifdef TUNE2
      if(lev.eq.1) then
        denlev=rho(1)
      else
        denlev=0.5*(rho(lev-1)+rho(lev))
      endif
#else
      if(lev.eq.1) denlev=rho(1)
      if(lev.ne.1) denlev=0.5*(rho(lev-1)+rho(lev))
#endif
      wpthp(lev)=fs(lev)/(denlev*cpa)		!turbulent flux of th at level
      wpqp(lev)=fq(lev)/(denlev*l)		!turbulent flux of q at level
 31   continue
      zoverh=z(lbase)/z(lbase+1)				!here assume lcl is base of cloud layer
      pmidmx=0.5*(p(lbase+1)+p(1))				!pressure at middle of mixed layer
      thovert=(pref/pmidmx)**gam				!theta/temp for mixed layer
#ifdef TUNE2
      tt1=thovert*l/cpa
      wpthep0=wpthp(1)+tt1*wpqp(1)			!turbulent flux of theta_e at surface
      wpthept=wpthp(lbase+1)+tt1*wpqp(lbase+1)		!turbulent flux of theta_e at top of b.l. cloud layer
#else
      wpthep0=wpthp(1)+thovert*l*wpqp(1)/cpa			!turbulent flux of theta_e at surface
      wpthept=wpthp(lbase+1)+thovert*l*wpqp(lbase+1)/cpa	!turbulent flux of theta_e at top of b.l. cloud layer
#endif

!	3.2) now calculate the solar flux terms from deardorff, 1976

      sumabs=0.0
      sumrho=0.0
      do 32 lay=1,lbase
      sumabs=sumabs+solin(lay)				!sum solar absorption in mixed layer (W/m**2)
      sumrho=sumrho+rho(lay)*dp(lay)
 32   continue
      rhomx=sumrho/(p(1)-p(lbase+1))			!mean density in mixed layer
#ifdef TUNE2
      tt1=1./(rhomx*cpa)
      rshmrso=-sumabs*tt1
      rscmrso=-(sumabs-solin(lbase))*tt1
#else
      rshmrso=-sumabs/(rhomx*cpa)
      rscmrso=-(sumabs-solin(lbase))/(rhomx*cpa)
#endif

!	3.3) calculate the turbulent buoyancy flux below (wpthvpb) and
!            above (wpthvpa) cloud base. then get mean buoyancy flux
!            in mixed layer (wpthvpm)

      de1=1.005							!values from Deardorff, 1976
      de2=8.3							!value from Deardorff, 1976
      de3=0.05							!value from Deardorff, 1976
      de4=0.5							!value from Deardorff, 1976
      wpqwp0=0.0						!assume no flux of cloud water
      wpqwpt=0.0						!assume no flux of cloud water
#ifdef TUNE2
      tt1=(1.-zoverh)*wpthep0 + zoverh*wpthept +
     1            zoverh*thovert*rshmrso - thovert*rscmrso 
      tt2=thsf*((1.-zoverh)*wpqwp0 + zoverh*wpqwpt)
      wpthvpb=de1*tt1 - de2*tt2
      wpthvpa=de4*tt1 - tt2
#else
      wpthvpb=de1*((1.-zoverh)*wpthep0 + zoverh*wpthept + 
     1            zoverh*thovert*rshmrso - thovert*rscmrso) -
     2            de2*thsf*((1.-zoverh)*wpqwp0 + zoverh*wpqwpt)

      wpthvpa=de4*((1.-zoverh)*wpthep0 + zoverh*wpthept +
     1           zoverh*thovert*rshmrso - thovert*rscmrso) -
     2           thsf*((1.-zoverh)*wpqwp0 + zoverh*wpqwpt) 
#endif

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
#ifdef TUNE2
!     p<4000Pa(=40hPa), cloud water may not be there...
      if (p1.lt.4000.0) then
        dqsdpma=0.0
      else
#endif
      p2=pk(i)+1000.0
      t1=t(i)
      th1=t1*(pref/p1)**gam
      call getthc2s(p1   ,p2   ,th1  ,th2  ,l    ,
     1              cpa  ,gam  ,pref )
      t2=th2*(p2/pref)**gam
      qs1=qsat(t1,p1)
      qs2=qsat(t2,p2)
      dqsdpma=(qs2-qs1)/(p2-p1)
#ifdef TUNE2
      endif
#endif
      p3=pk(i)
      t3=t(i)+0.5
      qs3=qsat(t3,p3)
      dqsdtp=(qs3-qs1)/(t3-t1)

!	4.2) calculate temperature tendency due to diabatic 
!	     processes (dtdtd)

      dtdtd = hrtr(i) + hrsl(i) + dtsur(i) + dt_ent(i) 

!	4.3) calculate saturation humidity tendency (dqsdt)

#ifndef TUNE1
      dqsdt=dqsdpma*(wpa(i)+0.05*g*cmf(i)) + dqsdtp*dtdtd	!assume 0.05 ratio between cmf and
#else
      dqsdt=dqsdpma*(wpa(i)+     g*cmf(i)) + dqsdtp*dtdtd	!assume 0.05 ratio between cmf and
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
#ifdef TUNE2
        tt1=lw(i)/delt
        tt2=tt1/e(i)
        e1 = e1 * tt2
        e2 = e2 * tt2
        e(i) = tt1
#else
        e1 = e1 * lw(i)/delt/e(i)
        e2 = e2 * lw(i)/delt/e(i)
        e(i) = lw(i)/delt
#endif
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
#ifdef TUNE2
        tt1=e2/lw(i)
        disp(i)=tt1*a(i)
        bb1(i)=tt1
#else
        disp(i)=e2*a(i)/lw(i)
        bb1(i)=e2/lw(i)			!save for analytic integration
#endif
      else				
        disp(i)=0.0
        bb1(i)=0.0			!save for analytic integration
      endif

 40   continue

      call precip(aa1,lw,a,t,q,rho,pk,dp,delt,lf,g,cpa,
     1            prate,lprate,fprate,gp,gpl,gpi,dd1,ep,
     2            epl,epi,dtmelt,dtprec,fliqex,ficeex,iflag)

!	update cloud fraction, liquid water/ice, temperature and humidity

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

!	partition source terms using prescribed fliq function
!	partition sink terms using existing liquid and ice ratios

      fl = fliq(t(i))
#ifdef TUNE2
      tt1=(1.0_r8 - fl)
      scv_l = fl*scv(i)
      scv_i = tt1*scv(i)
      sbl_l = fl*sbl(i)
      sbl_i = tt1*sbl(i)
      cond_l = fl*cond(i)
      cond_i = tt1*cond(i)
      e_l = fliqex(i)*e(i)
      e_i = ficeex(i)*e(i)
      if(dl_ent(i) .lt. 0.0) then			! here dl_ent is a sink
        dlent_l = fliqex(i)*dl_ent(i)
        dlent_i = ficeex(i)*dl_ent(i)
      else						! here dl_ent is a source
        dlent_l = fl*dl_ent(i)
        dlent_i = tt1*dl_ent(i)
      endif
#else
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
#endif

!	this first section uses finite differences

!	   cloud fraction:
!		- update cloud fraction a(i)
!		- keep a(i) between 0 and 1

!	   liquid water amount:
!		- update liquid water lw(i)
!		- keep lw(i) positive; avoid small values of lw(i)
!		- if lw(i)=0, set a(i) also to 0

      clwb(i) = clw(i)
      ciwb(i) = ciw(i)

      ifd=0							!ifd=0 ==> use analytic integration
      if(ifd.eq.0) goto 501

      a(i)=a(i)+delt*(sacv(i)+sabl(i)+sac(i)-disp(i))		!update a(i) using finite differences
      a(i)=min(1.0,a(i))
      a(i)=max(0.0,a(i))
      if(a(i).lt.1.0e-05) a(i)=0.0

!      lw(i)=lw(i)+delt*(scv(i)+sbl(i)+cond(i)-e(i)-gp(i)+dl_ent(i))	!finite differences

!     if(iflag.eq.1.and.i.eq.8) then
!       write(6,*)'starting clw(8)=',clw(8)
!       write(6,*)'starting ciw(8)=',ciw(8)
!       write(6,*)'starting lw(8)=',lw(8)
!       write(6,*)'dt*gp=',delt*gp(8)
!       write(6,*)'fliq+fice=',fliqex(8)+ficeex(8)
!     endif

!     if(iflag.eq.1) then
!       sumql0 = sumql0 + clw(i)*dp(i)/g
!       sumqi0 = sumqi0 + ciw(i)*dp(i)/g
!     endif

!     if(iflag.eq.1.and.i.eq.5) then
!       write(6,*)'before updates to clw, ciw:'
!       write(6,*)'  clw(5)=',clw(5)
!       write(6,*)'  ciw(5)=',ciw(5)
!       write(6,*)'    e_l =',e_l
!       write(6,*)'    e_i =',e_i
!       write(6,*)'  scv_l =',scv_l
!       write(6,*)'  scv_i =',scv_i
!       write(6,*)'  sbl_l =',sbl_l
!       write(6,*)'  sbl_i =',sbl_i
!       write(6,*)'dlent_l =',dlent_l
!       write(6,*)'dlent_i =',dlent_i
!     endif

!     clw(i) = clw(i) + delt*(scv_l + sbl_l + cond_l - e_l - gpl(i) + dlent_l)
!     ciw(i) = ciw(i) + delt*(scv_i + sbl_i + cond_i - e_i - gpi(i) + dlent_i)
      clw(i) = clw(i) + delt*(scv_l + sbl_l + cond_l - e_l - 
     1                        fliqex(i)*gp(i) + dlent_l)
      if(clw(i).lt.0.0) then		!assumes this happens when both e_l and gp*fliqex exceed supply
!       if(iflag.eq.1) write(6,*)'clw lt 0 at i=',i
        e_l = e_l + clw(i)/delt
        clw(i) = 0.0
      endif
      ciw(i) = ciw(i) + delt*(scv_i + sbl_i + cond_i - e_i - 
     1                        ficeex(i)*gp(i) + dlent_i)
      if(ciw(i).lt.0.0) then		!assumes this happens when both e_i and gp*ficeex exceed supply
!       if(iflag.eq.1) write(6,*)'ciw lt 0 at i=',i
        e_i = e_i + ciw(i)/delt
        ciw(i) = 0.0
      endif

!     if(iflag.eq.1) then
!       sumtest = sumtest + scv_i*delt*dp(i)*lf/g
!     endif

      lw(i) = clw(i) + ciw(i)

!     if(lw(i).lt.1.0e-10) then
!       dqfix=lw(i)
!       lw(i)=0.0
!     else
!       dqfix=0.0
!     endif
      if(lw(i).eq.0.0) a(i)=0.0
      goto 502

!	this section uses analytic integration - tiedtke, 1993
!
!	   cloud fraction:
!		- sum terms for eq(37) tiedtke 1993
!		- calculate new cloud fraction a(i)
!		- keep a(i) between 0 and 1
!
!	   liquid water amount:
!		- sum terms for eq(39) tiedtke 1993
!		- update liquid water lw(i)
!		- add effects of entrainment dl_ent(i)
!		- keep lw(i) positive; avoid small values of lw(i)
!		- if lw(i)=0, set a(i) also to 0
!		- calculate average conversion rate (gpnew) (eq(40) tiedtke, 1993)
!		- keep running sum for precipitation rate (pratenew)

 501  aa=aa1(i)+aa2(i)+aa3(i)
      bb=bb1(i)
      if(aa+bb.eq.0.0) then
        a(i)=a(i) + aa*delt
      else if((aa+bb)*delt.gt.50.0) then
        a(i)=aa/(aa+bb)
      else
#ifdef TUNE2
        tt1=exp(-(aa+bb)*delt)
        a(i)=a(i)*tt1 + (1.0-tt1)*aa/(aa+bb)	!update a(i) using analytic integration
#else
        a(i)=a(i)*exp(-(aa+bb)*delt) + 		!update a(i) using analytic integration
     1                  (1.0-exp(-(aa+bb)*delt))*aa/(aa+bb) 	
#endif
      endif
      a(i)=min(1.0,a(i))
      a(i)=max(0.0,a(i))
      if(a(i).lt.1.0e-05) a(i)=0.0

!     cc=scv(i)+sbl(i)+cond(i)-e(i)
!     dd=dd1(i)

      cc_l = scv_l + sbl_l + cond_l - e_l
      cc_i = scv_i + sbl_i + cond_i - e_i
      cc = cc_l + cc_i

      if(clw(i).eq.0.0) then
        dd_l = 0.0
      else
        dd_l = fliqex(i)*gp(i)/clw(i)
      endif

      if(ciw(i).eq.0.0) then
        dd_i = 0.0
      else
        dd_i = ficeex(i)*gp(i)/ciw(i)
      endif
      dd = dd_l + dd_i

      oldlw=lw(i)

!     if(iflag.eq.1.and.i.eq.2) then
!       write(6,*)'i = ',i
!       write(6,*)'clw, before =',clw(i)
!       write(6,*)'ciw, before =',ciw(i)
!       write(6,*)'cc_l=',cc_l
!       write(6,*)'  scv_l=',scv_l
!       write(6,*)'  sbl_l=',sbl_l
!       write(6,*)'  cond_l=',cond_l
!       write(6,*)'  e_l=',e_l
!       write(6,*)'cc_i=',cc_i
!       write(6,*)'  scv_i=',scv_i
!       write(6,*)'  sbl_i=',sbl_i
!       write(6,*)'  cond_i=',cond_i
!       write(6,*)'  e_i=',e_i
!       write(6,*)'dd_l=',dd_l
!       write(6,*)'dd_i=',dd_i
!       write(6,*)'dlent_l=',dlent_l
!       write(6,*)'dlent_i=',dlent_i
!       write(6,*)'gp=',gp(i)
!       write(6,*)'fliqex=',fliqex(i)
!       write(6,*)'epl=',epl(i)
!       write(6,*)'lpratenew=',lpratenew
!       write(6,*)'fpratenew=',fpratenew
!       write(6,*)'t=',t(i)
!       write(6,*)'dp=',dp(i)
!     endif


      if(dd_l .eq. 0.0) then
        clw(i) = clw(i) + cc_l*delt
      else if(dd_l*delt.gt.50.0) then
        clw(i) = cc_l/dd_l
      else
#ifdef TUNE2
        tt1=exp(-dd_l*delt)
        clw(i) = clw(i)*tt1 + cc_l*(1.0-tt1)/dd_l	! analytic integration
#else
        clw(i) = clw(i)*exp(-dd_l*delt) + 
     1           cc_l*(1.0-exp(-dd_l*delt))/dd_l	! analytic integration
#endif
      endif
      clw(i) = clw(i) + delt*dlent_l

      if(clw(i).lt.0.0) then		! assumes this happens when e_l exceeds supply
        e_l = e_l + clw(i)/delt		! reduce evaporation
        clw(i) = 0.0			! set clw=0
      endif

      if(dd_i .eq. 0.0) then
        ciw(i) = ciw(i) + cc_i*delt
      else if(dd_i*delt.gt.50.0) then
        ciw(i) = cc_i/dd_i
      else
#ifdef TUNE2
        tt1=exp(-dd_i*delt)
        ciw(i) = ciw(i)*tt1 + cc_i*(1.0-tt1)/dd_i	! analytic integration
#else
        ciw(i) = ciw(i)*exp(-dd_i*delt) + 
     1           cc_i*(1.0-exp(-dd_i*delt))/dd_i	! analytic integration
#endif
      endif
      ciw(i) = ciw(i) + delt*dlent_i

      if(ciw(i).lt.0.0) then		! assumes this happens when e_i exceeds supply
        e_i = e_i + ciw(i)/delt		! reduce evaporation
        ciw(i) = 0.0			! set ciw=0
      endif

      lw(i) = clw(i) + ciw(i)

      if(lw(i).eq.0.0) a(i)=0.0

!	need to recalculate gp and prate (lprate + fprate)

      gp(i) = (oldlw-lw(i))/delt + cc					! replace gp(i) with average net conversion rate
#ifdef TUNE2
      if(gp(i) .lt. 0.0.or.dd .eq. 0.0) gp(i)=0.0
#else
      if(gp(i) .lt. 0.0) gp(i)=0.0
      if(dd .eq. 0.0) gp(i)=0.0
#endif
      if(t(i) .le. 273.15) then					! here all ice
        gpl(i) = 0.0						! no liquid precip
        gpi(i) = gp(i)						! all precip is ice
        dtprec(i) = fliqex(i)*gp(i)*lf*delt/cpa			! if some frozen precip formed from liquid, account for dT
      else							! else, all liquid
        gpl(i) = gp(i)						! here all precip is liquid
        gpi(i) = 0.0						! no frozen precip
        dtprec(i) = -ficeex(i)*gp(i)*lf*delt/cpa		! if some liquid precip formed from ice, account for dT
      endif

#ifdef TUNE2
cc      tt1=dp(i)/g
      tt1=dpg(i)
      if(epl(i)*tt1 .gt. lpratenew) epl(i)=lpratenew/tt1	! make sure evap. does not exceed available lpratenew
      if(epi(i)*tt1 .gt. fpratenew) epi(i)=fpratenew/tt1	! make sure evap. does not exceed available fpratenew

      lpratenew = lpratenew + (gpl(i)-epl(i))*tt1
      lpratenew = max(0.0,lpratenew)
      fpratenew = fpratenew + (gpi(i)-epi(i))*tt1
      fpratenew = max(0.0,fpratenew)
#else
      if(epl(i)*dp(i)/g .gt. lpratenew) epl(i)=lpratenew/(dp(i)/g)	! make sure evap. does not exceed available lpratenew
      if(epi(i)*dp(i)/g .gt. fpratenew) epi(i)=fpratenew/(dp(i)/g)	! make sure evap. does not exceed available fpratenew

      lpratenew = lpratenew + (gpl(i)-epl(i))*dp(i)/g
      lpratenew = max(0.0,lpratenew)
      fpratenew = fpratenew + (gpi(i)-epi(i))*dp(i)/g
      fpratenew = max(0.0,fpratenew)
#endif

      pratenew = lpratenew + fpratenew

!     pratenew=pratenew+(gpnew-ep(i))*(-dz(i))*rho(i)
!     if(pratenew.lt.0.0) then
!       dqfix=dqfix+pratenew*delt/(-dz(i)*rho(i))
!       pratenew=0.0
!     endif

!	now calculate changes in the humidity and temperature

!502  dqliw(i)=dqfix+delt*(-sbl(i)-cond(i)+e(i)+ep(i)+dq_ent(i))
#ifdef TUNE2
 502  dqliw(i) = dqfix + delt*
     *                       ((-sbl_l - cond_l + e_l + epl(i)) + 
     1                        (-sbl_i - cond_i + e_i + epi(i)) + 
     2                        dq_ent(i) )
#else
 502  dqliw(i) = dqfix + delt*(-sbl_l - cond_l + e_l + epl(i)) + 
     1                   delt*(-sbl_i - cond_i + e_i + epi(i)) + 
     2                   delt*dq_ent(i)
#endif

!     if(iflag.eq.1.and.i.eq.2) then
!       write(6,*)'layer=',i
!       write(6,*)'q=',q(i)
!       write(6,*)'sbl_l,sbl_i=',delt*sbl_l,delt*sbl_i
!       write(6,*)'cond_l,cond_i=',delt*cond_l,delt*cond_i
!       write(6,*)'e_l,e_i=',delt*e_l,delt*e_i
!       write(6,*)'epl,epi=',delt*epl(i),delt*epi(i)
!       write(6,*)'dq_ent(i)=',delt*dq_ent(i)
!     endif

!     if(iflag.eq.1.and.i.eq.2) then
!     if(iflag.eq.1) then
!       write(6,*)'  '
!       write(6,3993)i,clw(i), delt*cond_l, delt*e_l, delt*fliqex(i)*gp(i)
!       write(6,3994)i,ciw(i), delt*cond_i, delt*e_i, delt*ficeex(i)*gp(i)
!       write(6,3995)i,q(i),delt*(cond_l+cond_i),delt*(e_l+e_i),delt*(epl(i)+epi(i)),q(i)+dqliw(i)
!3993   format('k= ',i3,'  clw=',1pe10.2,'  dt*cond_l=',e10.2,'  dt*e_l=',e20.12,'  dt*fliq*gp=',e10.2)
!3994   format('k= ',i3,'  ciw=',1pe10.2,'  dt*cond_i=',e10.2,'  dt*e_i=',e20.12,'  dt*fice*gp=',e10.2)
!3995   format('k= ',i3,'    q=',1pe10.2,'  dt*cond  =',e10.2,'  dt*e  =',e20.12,'  dt*epl    =',e10.2,'  dqnew=',e10.2)
 
!       sumgpr = sumgpr + delt*(fliqex(i)+ficeex(i))*gp(i)*dp(i)/g
!       sumepr = sumepr + delt*(epl(i)+epi(i))*dp(i)/g

!       sumqv0 = sumqv0 + q(i)*dp(i)/g
!       sumqv1 = sumqv1 + (q(i)+dqliw(i))*dp(i)/g
!       sumql1 = sumql1 + clw(i)*dp(i)/g
!       sumqi1 = sumqi1 + ciw(i)*dp(i)/g
!       sumfix = sumfix + dqfix*dp(i)/g

!       sumqtot = 0.0
!       do k = 1,nl
!         if(k.ge.i) then
!           sumqtot = sumqtot + (q(k)+dqliw(k)+clw(k)+ciw(k))*dp(k)/g		! already updated, include dqliw but not dlf
!         else
!           sumqtot = sumqtot + (q(k)+clw(k)+ciw(k)+dlf(k)*delt)*dp(k)/g	! not updated, include dlf, but not dqliw
!         endif
!       enddo
!       sumqtot = sumqtot + prate*delt
!       write(6,*)'after layer =',i,'  sumqtot=',sumqtot

!     endif

#ifdef TUNE2
      dtliw(i) = delt*
     *              ( ((sbl_l + cond_l - e_l - epl(i))*l/cpa) + 
     1                ((sbl_i + cond_i - e_i - epi(i))*(l+lf)/cpa) + 
     2                dt_ent(i) + 
     3                scv_i*lf/cpa)					! this last term is for detrained liquid water freezing
#else
      dtliw(i) = delt*((sbl_l + cond_l - e_l - epl(i))*l/cpa) + 
     1           delt*((sbl_i + cond_i - e_i - epi(i))*(l+lf)/cpa) + 
     2           delt*dt_ent(i) + 
     3           delt*scv_i*lf/cpa					! this last term is for detrained liquid water freezing
#endif

      dtliw(i) = dtliw(i) + dtmelt(i) + dtprec(i)

!     if(t(i).ge.273.15) then
!      dtliw(i)=delt*(l*(sbl(i)+cond(i)-e(i)-ep(i))/cpa+dt_ent(i))
!      if(iflag.eq.1) then
!         sume=sume-delt*l*e(i)*dp(i)/g
!         sumep=sumep-delt*l*ep(i)*dp(i)/g
!         sumcnd=sumcnd+delt*l*cond(i)*dp(i)/g
!      endif
!     else
!      dtliw(i)=delt*((l+lf)*(sbl(i)+cond(i)-e(i)-ep(i))/cpa+dt_ent(i))
!      if(iflag.eq.1) then
!         sume=sume-delt*(l+lf)*e(i)*dp(i)/g
!         sumep=sumep-delt*(l+lf)*ep(i)*dp(i)/g
!         sumcnd=sumcnd+delt*(l+lf)*cond(i)*dp(i)/g
!      endif
!     endif
!     dtliw(i)=dtliw(i)+dtmelt(i)
!     if(iflag.eq.1) then
!       summlt = summlt + dtmelt(i)*cpa*dp(i)/g
!       sumfix = sumfix + dqfix*dp(i)/g
!       sumtot = sumtot + dtliw(i)*cpa*dp(i)/g
!     endif

 500  continue
!     if(iflag.eq.1) then
!       write(6,*)'             sumtest=',sumtest,'  units = J/m**2)'
!       write(6,*)'             sume  =',sume
!       write(6,*)'             sumep =',sumep
!       write(6,*)'             sumcnd=',sumcnd
!       write(6,*)'             summlt=',summlt
!       write(6,*)'             sumfix=',sumfix
!       write(6,*)'             sumtot=',sumtot
!       write(6,*)'             sumepr=',sumepr
!       write(6,*)'             sumgpr=',sumgpr
!       write(6,*)'             sumqv0=',sumqv0
!       write(6,*)'             sumql0=',sumql0
!       write(6,*)'             sumqi0=',sumqi0
!       write(6,*)'             total0=',sumqv0+sumql0+sumqi0
!       write(6,*)'             sumqv1=',sumqv1
!       write(6,*)'             sumql1=',sumql1
!       write(6,*)'             sumqi1=',sumqi1
!       write(6,*)'             total1=',sumqv1+sumql1+sumqi1
!     endif

      if(ifd.eq.0) then
        prate=pratenew
        lprate=lpratenew
        fprate=fpratenew
      endif

      wvcol=0.0
      wlcol=0.0
      wicol=0.0
      dqcol=0.0
      dqcolabs=0.0
      do 68 i=1,nl
#ifdef TUNE2
        tt1=dpg(i)
        dqcol = dqcol + dqliw(i)*tt1
        dqcolabs = dqcolabs + abs(dqliw(i))*tt1
        wvcol = wvcol + q(i)*tt1
        wlcol = wlcol + clw(i)*tt1
        wicol = wicol + ciw(i)*tt1
#else
        dqcol = dqcol + dqliw(i)*dp(i)/g
        dqcolabs = dqcolabs + abs(dqliw(i))*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
#endif
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
#ifdef TUNE2
        tt1=dpg(i)
        dqbef = dqliw(i)
!       dqliw(i) = (1.0_r8 - dqadjr) * dqliw(i)
        if(dqcolabs .ne. 0.0) then
          dqliw(i) = dqliw(i) - 
     1               (dqadjr*abs(dqliw(i)*tt1)/dqcolabs)/tt1
!       else
!         dqliw(i) = dqliw(i) - dqadjr*(q(i)*tt1)/qcol
        endif
        dqcol = dqcol + dqliw(i)*tt1
        wvcol = wvcol + q(i)*tt1
        wlcol = wlcol + clw(i)*tt1
        wicol = wicol + ciw(i)*tt1
#else
        dqbef = dqliw(i)
!       dqliw(i) = (1.0_r8 - dqadjr) * dqliw(i)
        if(dqcolabs .ne. 0.0) then
          dqliw(i) = dqliw(i) - 
     1               (dqadjr*abs(dqliw(i)*dp(i)/g)/dqcolabs)*g/dp(i)
!       else
!         dqliw(i) = dqliw(i) - dqadjr*(q(i)*dp(i)/g)/qcol
        endif
        dqcol = dqcol + dqliw(i)*dp(i)/g
        wvcol = wvcol + q(i)*dp(i)/g
        wlcol = wlcol + clw(i)*dp(i)/g
        wicol = wicol + ciw(i)*dp(i)/g
#endif
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
#ifdef TUNE2
        tt1=dpg(i)
        tt2=cpa*tt1
        dtcol = dtcol + dtliw(i)*tt2
        dtcolabs = dtcolabs + abs(dtliw(i))*tt2
        dqcol = dqcol + dqliw(i)*tt1
        secol = secol + cpa*t(i)*tt1
#else
        dtcol = dtcol + dtliw(i)*cpa*dp(i)/g
        dtcolabs = dtcolabs + abs(dtliw(i))*cpa*dp(i)/g
        dqcol = dqcol + dqliw(i)*dp(i)/g
        secol = secol + cpa*t(i)*dp(i)/g
#endif
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
#ifdef TUNE2
        tt1=cpa*dpg(i)/g
        if(dtcolabs .ne. 0.0) then
          dtliw(i) = dtliw(i) - 
     1    (dtadjr*abs(dtliw(i)*tt1)/dtcolabs)/tt1
        endif
        dtcol = dtcol + dtliw(i)*tt1
#else
!       dtliw(i) = (1.0_r8 - dtadjr) * dtliw(i)
        if(dtcolabs .ne. 0.0) then
          dtliw(i) = dtliw(i) - 
     1    (dtadjr*abs(dtliw(i)*cpa*dp(i)/g)/dtcolabs)*g/dp(i)/cpa
        endif
        dtcol = dtcol + dtliw(i)*cpa*dp(i)/g
#endif
 79   continue
      te_end = dtcol + secol +(l+lf)*(wvcol+dqcol) + lf*wlcol + 
     1         lprate*delt*lf
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

#ifdef TUNE2
      real(r8) :: tt1
      real(r8) :: tt2
#endif

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
#ifdef TUNE2
      tt1=(pmidc/pref)**gam
      tt2=thovertmid*l/cpa
      do it=1,100
      tcg=thcg*tt1
      qsg=qsat(tcg,pmidc)
      tcmid=0.5*(tcbase+tcg)
      qsmid=0.5*(qsbase+qsg)
      dqs=qsg-qsbase
      dtc=tcg-tcbase
      thcnew=thcbase-tt2*(dqs-(qsmid*dtc/tcmid))
      if(abs(thcnew-thcg).gt.0.01) then
        thcg=0.25*thcnew+0.75*thcg
      else
        exit
      endif
      end do
#else
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
#endif
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

#ifdef TUNE2
! tt1,tt2,phalf,thovertmid are already calculated
      thcg=thctop-1.0			!initial guess for cloud potential temp.
      do it=1,100
      tcg=thcg*tt1       
      qsg=qsat(tcg,pmidc)
      tcmid=0.5*(tctop+tcg)
      qsmid=0.5*(qstop+qsg)
      dqs=qsg-qstop
      dtc=tcg-tctop
      thcnew=thctop-tt2*(dqs-(qsmid*dtc/tcmid))
      if(abs(thcnew-thcg).gt.0.01) then
        thcg=0.25*thcnew+0.75*thcg
      else
        exit
      endif
      end do
#else
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
#endif
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
#ifdef TUNE2
      real(r8) :: tt1
      real(r8) :: tt2
#endif

      integer :: it

      it=0				!initialize iteration counter
      phalf=0.5*(p2+p1)			!pressure halfway between cloud base and mid cloud height
#ifdef TUNE2
cc      thcg=th1-0.5			!initial guess for cloud potential temp.
      if (p1.gt.30000.0) then
        thcg=th1-0.5			!initial guess for cloud potential temp.
      else
        thcg=th1-0.1			!initial guess for cloud potential temp.
      endif
#else
      thcg=th1+1.0			!initial guess for cloud potential temp.
#endif
      thovertmid=((pref/phalf)**gam)
      t1=th1*(p1/pref)**gam
      qs1=qsat(t1,p1)

#ifdef TUNE2
      tt1=(p2/pref)**gam
      tt2=thovertmid*l/cpa
      do it=1,100
      tcg=thcg*tt1
      qsg=qsat(tcg,p2)
      tcmid=0.5*(th1+tcg)
      qsmid=0.5*(qs1+qsg)
      dqs=qsg-qs1
      dtc=tcg-t1
      thcnew=th1-tt2*(dqs-(qsmid*dtc/tcmid))
      if(abs(thcnew-thcg).gt.0.01) then
        thcg=0.25*thcnew+0.75*thcg
      else
        exit
      endif
      end do
#else
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
#endif
      th2=0.5*(thcnew+thcg)		!cloud potential temperature at p=pmidc
c****
c      write(*,*) '****'
c      write(*,*) it
c      write(*,*) p1,p2
c      write(*,*) th1,th2
c****

      return
      end
!************************************************************************************
      subroutine precip(aa1   ,lw    ,a     ,t     ,q     ,
     1                  rho   ,pk    ,dp    ,delt  ,lf    ,
     2                  g     ,cpa   ,prate ,lprate,fprate,
     3                  gp    ,gpl   ,gpi   ,dd1   ,ep    ,
     4                  epl   ,epi   ,dtmelt,dtprec,fliqex,ficeex,iflag)

      use shr_kind_mod,  only: r8=>shr_kind_r8
      use ppgrid,        only: nl=>pver

      implicit none
      
!	input quantities

      real(r8), intent(in) :: aa1(nl)		! term containing convective detrainment
      real(r8), intent(in) :: lw(nl)		! grid-mean cloud water (liquid + ice) (kg/kg)
      real(r8), intent(in) :: a(nl)		! cloud fraction
      real(r8), intent(in) :: t(nl)		! temperature
      real(r8), intent(in) :: q(nl)		! specific humidity
      real(r8), intent(in) :: rho(nl)		! air density
      real(r8), intent(in) :: pk(nl)		! pressure at midlevels (Pa)
      real(r8), intent(in) :: dp(nl)		! layer pressure thickness (Pa)
      real(r8), intent(in) :: fliqex(nl)	! liquid fraction of existing cloud water
      real(r8), intent(in) :: ficeex(nl)	! ice fraction of existing cloud water
      real(r8), intent(in) :: delt		! model timestep (sec)
      real(r8), intent(in) :: lf		! latent heat of fusion
      real(r8), intent(in) :: g			! gravity
      real(r8), intent(in) :: cpa		! 

!	output quantities

      real, intent(out) :: prate		! precipitation rate at surface
      real, intent(out) :: fprate		! frozen precipitation rate at surface
      real, intent(out) :: lprate		! liquid precipitation rate at surface
      real, intent(out) :: gp(nl)		! conversion of cloud water to precipitation (kg/kg/sec)
      real, intent(out) :: gpl(nl)		! conversion of liquid cloud water to precipitation (kg/kg/sec)
      real, intent(out) :: gpi(nl)		! conversion of ice cloud water to precipitation (kg/kg/sec)
      real, intent(out) :: dd1(nl)		! conversion of cloud water to precipitation for analytic integration
      real, intent(out) :: ep(nl)		! evaporation of precipitation falling into layer (kg/kg/sec)
      real, intent(out) :: epl(nl)		! evaporation of liquid precipitation falling into layer (kg/kg/sec)
      real, intent(out) :: epi(nl)		! evaporation of frozen precipitation falling into layer (kg/kg/sec)
      real, intent(out) :: dtmelt(nl)		! temperature change due to melting of frozen precipitation
      real, intent(out) :: dtprec(nl)		! temperature change due to formation of liquid or ice precipitation

!	local derived quantities

      real(r8) :: ap				! fractional area where precipiation is occurring
      real(r8) :: c0				! adjusted characteristic conversion time
      real(r8) :: c0star			! characteristic conversion time
      real(r8) :: lcrit				! adjusted critical cloud water for precipitation
      real(r8) :: lcrits			! critical cloud water for precipitation 
      real(r8) :: lwc				! in-cloud cloud water content (liquid + ice)
      real(r8) :: f1				! coefficient for coalescense processes
      real(r8) :: f2				! coefficient for Bergeron-Findeisen processes
      real(r8) :: fci_l				! coefficient for cold-cirrus microphysical properties
      real(r8) :: fci_c				! coefficient for cold-cirrus microphysical properties
      real(r8) :: ss, x, y, prp			! temporary storage
      real(r8) :: fl				! liquid fraction of cloud
      real(r8) :: fliq				! real(r8) :: function
      real(r8) :: qsat				! real(r8) :: function
      real(r8) :: rhocw				! density of cloud water
      real(r8) :: rc				! characteristic droplet radius (Manton-Cotton)
      real(r8) :: vc				! characteristic terminal velocity (Manton-Cotton)
      real(r8) :: fc				! collision frequency (Manton-Cotton)
      real(r8) :: rhocwmin			! minimum cloud water density for coalescence to occur
      real(r8) :: sc				! conversion of cloud water to precip (Manton-Cotton)
      real(r8) :: del_a				! temporary storage
      real(r8) :: fac				! temporary storage
      real(r8) :: fprate_sav			! temporary storage 
      real(r8) :: dtmx				! temporary storage
      integer  :: i				! counter

      integer :: iflag

!	constants and parameters

      real, parameter :: c1 = 300.0		! constant used in derivation of coalescense coefficient (Sundqvist) (kg/(m**2 sec))**-0.5
      real, parameter :: c2 = 0.5		! constant used in derivation of Bergeron-Findeison coefficient (Sundqvist) (1/degK)
      real, parameter :: c0star_c = 1.0e-03	! characteristic conversion time: convective case (Sundqvist) (1/sec) 
      real, parameter :: c0star_s = 1.0e-04	! characteristic conversion time: stratiform case (Sundqvist) (1/sec)
      real, parameter :: lcrits_c = 5.0e-04	! critical cloud water for precipitation: convective case (Sundqvist) (kg/kg)
      real, parameter :: lcrits_s = 3.0e-04	! critical cloud water for precipitation: stratiform case (Sundqvist) (kg/kg)
      real, parameter :: dropcn = 50.0		! droplet concentration (Manton-Cotton) (1/cm**3) constant for now
      real, parameter :: rhowat = 1.0		! water density (Manton-Cotton) (g/cm**3)
      real, parameter :: visc = 1.72e-04	! viscosity of air (Manton-Cotton) (g/cm/sec)
      real, parameter :: ec = 0.55		! mean collector efficiency (Manton-Cotton)

#ifdef TUNE2
      real(r8) :: tt1
      real(r8) :: tt2
#endif

!	Initialize

      prate=0.0				!units = kg/m**2/sec
      lprate=0.0			!units = kg/m**2/sec
      fprate=0.0			!units = kg/m**2/sec
      ap=0.0				!fractional prec. area

      do 60 i=nl,1,-1

      if(aa1(i).gt.0.0) then		! if aa1(i) > then have convective detrainment of cloud water
        c0star=c0star_c
        lcrits=lcrits_c
      else
        c0star=c0star_s
        lcrits=lcrits_s
      endif

!	6.1) conversion of cloud water to prec (gp)

      if(lw(i).ne.0.0) then

        if(a(i).ne.0.0) then
          lwc=lw(i)/a(i)
        else
          lwc=0.0
        endif

!	coef for coalescence processes (f1) (fco in sundqvist)

        f1=1.0+c1*sqrt(prate)

!	coef for Bergeron-Findeisen processes (f2) (fbf in sundqvist)
#ifdef TUNE2
        if(t(i).ge.268.0) then
          f2=1.0
        else if(t(i).gt.250.0) then
          f2=1.0+c2*sqrt(268.0-t(i))
        else
          f2=1.0+c2*4.243                               !4.243=sqrt(268-250)
        endif
#else
        if(t(i).ge.268.0) f2=1.0
        if(t(i).lt.268.0.and.t(i).gt.250.0) f2=1.0+c2*sqrt(268.0-t(i))
        if(t(i).le.250.0) f2=1.0+c2*4.243				!4.243=sqrt(268-250)
#endif

!	coefficients for cold cirrus microphysical properties 
!	  fci_l - coef for lcrit
!	  fci_c - coef for c0

        if(t(i).le.250.0.and.aa1(i).eq.0.0) then		! aa1(i)=0 ==> stratiform cases
          x=abs(t(i)-233.0)/17.0
          y=x*(1.0+x+4.0*x*x/3.0)
#ifdef TUNE2
          if(t(i).lt.233.0) then
            prp=-1.0
          else
            prp=1.0
          endif
#else
          if(t(i).lt.233.0) prp=-1.0
          if(t(i).ge.233.0) prp=1.0
#endif

          fci_l=1.0/(0.15*(1.07+prp*y/(1.0+y)))	!coef for cirrus clouds
          fci_l=fci_l/f2			!acts to replace f2 when t<250

          fci_c=min( 1.0+4.0*(250.0-t(i))/18.0 , 5.0 )
          
        else
          fci_l=1.0
          fci_c=1.0
        endif

        c0=c0star*f1*f2*fci_c
        lcrit=lcrits/(f1*f2*fci_l)
        ss=max( -80.0 , -((lwc/lcrit)**2) )

!	gp(i)=precipitation from this layer (units=kg/kg/sec)

#ifdef TUNE2
        tt1=c0*(1.0-exp(ss))
        gp(i)=lw(i)*tt1
        dd1(i)=tt1                                      !save for analytic integration
#else
        gp(i)=lw(i)*c0*(1.0-exp(ss))
        dd1(i)=c0*(1.0-exp(ss))				!save for analytic integration
#endif

!	alternative formulation of precipitation from manton and cotton 1977
!	use only for liquid water region

        fl = fliq(t(i))					!liquid fraction of cloud water
        if(fl.eq.0.0) goto 3290				!skip manton/cotton if no liquid water present

        rhocw = lwc*rho(i)				!density of cloud water (kg/m**3)
        rhocw = 0.001 * rhocw				!put into units of g/cm**3

        rc = ((0.75/3.14159)*rhocw/(rhowat*dropcn))**(1.0/3.0)	!characteristic droplet radius (cm)
        vc = 2.0*rhowat*980.0*rc*rc/(9.0*visc)		!characteristic terminal velocity (cm/sec)
        fc = 3.14159*rc*rc*ec*vc*dropcn			!collision frequency (1/sec) 

        rhocwmin = 0.00419*dropcn			!minimum cloud water for coalescence to occur (g/m**3)
        rhocwmin = 1.0e-06*rhocwmin			!units = g/cm**3
        if(rhocw.lt.rhocwmin) then
          sc = 0.0
        else
          sc = fc * rhocw				!conversion of cloud water to rain (g/cm**3/sec)
        endif

        sc = 1000.0*sc					!units = kg/m**3/sec
        sc = sc/rho(i)					!units = kg_w/kg_a/sec
        sc = sc*a(i)					!multiply by cloud fraction to get grid-mean value

!	replace tiedtke/sundqvist autoconversion with manton-cotton

!       gp(i) = fl*sc + (1.0-fl)*lw(i)*c0*(1.0-exp(ss))
!       dd1(i) = gp(i)/lw(i)

 3290   if(gp(i).lt.0.0) then
          gp(i)=0.0
          dd1(i)=0.0
        endif

!	make sure gp(i) does not exceed available 
!	liquid water in this layer

        if(gp(i)*delt.gt.lw(i)) then
          gp(i)=lw(i)/delt
          dd1(i)=1.0/delt
        endif
        if(gp(i).lt.0.0) then
          gp(i)=0.0
          dd1(i)=0.0
        endif

!	precip. forming is either all ice or liquid depending on temp

        if(t(i) .le. 273.15) then			! here all ice
          gpl(i) = 0.0					! no liquid precip is formed
          gpi(i) = gp(i)				! all precip is ice
          dtprec(i) = fliqex(i)*gp(i)*lf*delt/cpa	! some frozen precip formed from liquid cloud water; account for dT
        else						! else all liquid
          gpl(i) = gp(i)				! all precip is liquid
          gpi(i) = 0.0					! no ice precip is formed
          dtprec(i) = -ficeex(i)*gp(i)*lf*delt/cpa	! some liquid precip formed from ice cloud water; account for dT
        endif

      else

        gp(i) = 0.0
        gpl(i) = 0.0
        gpi(i) = 0.0
        dd1(i) = 0.0
        dtprec(i) = 0.0

      endif

!	6.2) evaporation of rain falling into layer from above (ep)

      if(prate.eq.0.0) then
        ep(i) = 0.0
        epl(i) = 0.0
        epi(i) = 0.0
      else
        del_a=max(ap-a(i),0.0)
        if(ap.eq.0.0) then
          write(6,*)'PROBLEM WITH AP'
          write(6,*)'layer=',i
          write(6,*)'prate=',prate
!         stop
        endif
        fac=(sqrt(pk(i)/pk(1))*prate/(5.9e-03*ap))**0.577
        ep(i)=del_a*5.44e-04*(qsat(t(i),pk(i))-q(i))*fac
        epl(i) = ep(i)
        epi(i) = ep(i)
        if(epl(i)*dp(i)/g .gt. lprate) epl(i)=lprate/(dp(i)/g)
        if(epl(i).lt.1.0e-20) epl(i)=0.0
        if(epi(i)*dp(i)/g .gt. fprate) epi(i)=fprate/(dp(i)/g)
        if(epi(i).lt.1.0e-20) epi(i)=0.0
      endif

!	6.3) calculate net precipitation falling thru layer (prate)

!     prate=prate+(gp(i)-ep(i))*dp(i)/g
!     prate=max(0.0,prate)
      lprate = lprate + (gpl(i)-epl(i))*dp(i)/g
      lprate = max(0.0,lprate)
      fprate = fprate + (gpi(i)-epi(i))*dp(i)/g
      fprate = max(0.0,fprate)

      prate = lprate + fprate

      if(prate.gt.0.0) ap=max(ap,a(i))		!precipitation area above
 
!	6.4) keep track of whether precipitation is frozen or not

!     if(t(i).le.273.15) fprate=prate

!	6.4) check to see if any frozen precipitation melts or
!	     if any liquid precip freezes, and adjust dtmelt if necessary

      if(t(i).gt.273.15.and.fprate.gt.0.0) then
        dtmelt(i) = -fprate*lf*g*delt/cpa/dp(i)		! neg. value as ice melts
        dtmx = 273.15-t(i)				! maximum allowed temperature decrease (neg. value)
        if(dtmelt(i).lt.dtmx) then
          lprate = lprate + (dtmx/dtmelt(i))*fprate	! add some melted precipitation to liquid precip
          fprate = (1.0-dtmx/dtmelt(i))*fprate		! do not melt all precipitation
          dtmelt(i) = dtmx				! adjust temperature decrease due to melting
        else
          lprate = lprate + fprate			! add all melted precipitation to liquid precip
          fprate = 0.0					! all frozen precip has now melted
        endif
      endif

      if(t(i) .le. 273.15 .and. lprate .gt. 0.0) then
        dtmelt(i) = lprate*lf*g*delt/cpa/dp(i)		! pos value as liquid freezes
        dtmx = 273.15 - t(i)				! maximum allowed temperature increase
        if(dtmelt(i) .gt. dtmx) then
#ifdef TUNE2
          tt1=(dtmx/dtmelt(i))*lprate
          fprate = fprate + tt1				! add newly frozen
          lprate = lprate - tt1				! not all liquid precipitation freezes
#else
          fprate = fprate + (dtmx/dtmelt(i))*lprate	! add newly frozen 
          lprate = (1.0-dtmx/dtmelt(i))*lprate		! not all liquid precipitation freezes
#endif
          dtmelt(i) = dtmx				! adjust temperature increase due to freezing
        else
          fprate = fprate + lprate			! add all newly frozen precip to existing frozen precip
          lprate = 0.0					! all liquid precip has now frozen
        endif
      endif

      prate = lprate + fprate				! total should not change, but just to make sure

 60   continue

      return
      end

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
#ifdef TUNE2
      real(r8) :: tt1
      real(r8) :: tt2
#endif

#ifdef TUNE2
      tt1=1.0/(0.622+(0.378*q(lev)))
      tt2=log(6.11)
#endif

#ifdef TUNE2
      pg=pk(lev)-5000.0
#else
      pg=pk(lev)-20000.0
#endif
 300  do 400 n=1,500
#ifdef TUNE2
      pgmb=pg*1.e-2
      pgold=pg
      es=(pgmb*q(lev))*tt1
      lnterm=log(es)-tt2
#else
      pgmb=pg/100.0
      pgold=pg
      es=(pgmb*q(lev))/(0.622+(0.378*q(lev)))
      lnterm=log(es/6.11)
#endif
      tg=(4714.742-(35.86*lnterm))/(17.26-lnterm)
      pg=((tg/th(lev))**3.497)*pref
      diffpg=pg-pgold
#ifdef TUNE2
      if(abs(diffpg).lt.50.0) exit
#else
      if(abs(diffpg).lt.50.0) goto 1111
#endif
 400  continue
 1111 pval=pg
c****
cc      write(*,*) '****'
cc      write(*,*) n
cc      write(*,*) pk(lev),pg
c****
      if(pval.gt.pk(lev)) pval=pk(lev)
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
#ifdef TUNE2
      if(tc.ge.0.0) then
        fliq_result=1.0
      else if(tc.le.-15.0) then
        fliq_result=0.0
      else if(tc.ge.-5.0) then
        fliq_result=1.0-((0.2*tc)**2)/3.0
      else 
        fliq_result=(((tc+15.0)/5.0)**2.0)/6.0
      endif
#else
      if(tc.ge.0.0) fliq_result=1.0
      if(tc.lt.0.0.and.tc.ge.-5.0) fliq_result=1.0-((0.2*tc)**2)/3.0
      if(tc.lt.-5.0.and.tc.gt.-15.0) 
     1   fliq_result=(((tc+15.0)/5.0)**2.0)/6.0
      if(tc.le.-15.0) fliq_result=0.0
#endif

      return
      end

