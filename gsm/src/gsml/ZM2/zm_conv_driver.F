      subroutine zm_conv_driver(imx2,imx22,kmx,
     &              jcap,delt,si,ps,nstep,
     &              q1,t1,rn,kbot,ktop,kuo,lat,
     &              hpbl,gamt,gamq,qci,qrs,delx,
     &              tp, qp,
     &              tp1, qp1, clw, detr)
!   
      use ppgrid,        only: pcols, pver, pverp
      use zm_conv,       only: zm_convi, zm_convr, zm_conv_evap
      use physconst,     only: rair, gravit, cpair, zvir, latvap
      use physics_types, only: physics_state, physics_ptend
!   
c-----------------------------------------------------------------------
c   driver subroutine for 
c   modified Zhang McFarlane deep convection scheme
c
c   jul.2006  A.Shimpo
c     * created from ccmcnv_driver.F
c     * detrained cloud water evaporates for non cloud water scheme
c     * "define EVAP" : evaporation of convective precipitation
c       (default)
c
c-----------------------------------------------------------------------
#include <define.h>
#include <implicit.h>

#define EVAP  ! includes evaporation of convective precipitation

c------------------------------arguments--------------------------------
c
c input arguments
c
      integer nstep                 ! current time step for cmfmca diag.
      integer lat                   ! latitude index (s->n)
      real    ps(imx2)
      real    si(imx22,kmx+1)
      real    t1(imx22,kmx)
      real    q1(imx22,kmx)
      real    hpbl(imx2)
      real    gamt(imx2),gamq(imx2)
c
c input/output arguments
c
      real    tp(imx22,kmx)
      real    qp(imx22,kmx)
      real    tp1(imx22,kmx)
      real    qp1(imx22,kmx)
c
c output arguments
c
      integer kbot(imx2),ktop(imx2),kuo(imx2)
      real    qci(imx2,kmx),qrs(imx2,kmx)
      real    detr(imx22,kmx)         !
      real    clw(imx22,kmx)          ! [kg/kg]
      real    rn(imx2)                ! [m]

c
c local
c
      type(physics_state) :: state
      type(physics_ptend) :: ptend

      real precc(pcols)             ! convective-scale preciptn rate [m/s]
      real snow(pcols)              ! convective-scale snowfall rate [m/s]
      real cnt(pcols)               ! top level of convective activity
      real cnb(pcols)               ! bottom level of convective activity
      real cmfdqr(pcols,pver)       ! dq/dt due to moist convective rainout 
      real cmfmc(pcols,pverp)       ! moist convection cloud mass flux
      real psfc(pcols)              ! surface temperature
      real tv(pcols,pver)           ! virtual temperature array
      real dlf(pcols,pver)          ! detraining cld h20 from convection
      real pflx(pcols,pverp)        ! conv rain flux thru out btm of lev
      real pmln(pcols,pver)         ! log(pmid)
      real piln(pcols,pverp)        ! log(pint)
      real cldfrc(pcols,pver)       ! cloud fraction (evap)

      real zdu(pcols,pver)
      real rprd(pcols,pver)
      real mu(pcols,pver)
      real md(pcols,pver)
      real du(pcols,pver)
      real eu(pcols,pver)
      real ed(pcols,pver)
      real dp(pcols,pver)
      real dsubcld(pcols)
      real ql(pcols,pver)      
      real rliq(pcols)        
      real ptendq(pcols,pver)  
      real ptendt(pcols,pver) 
      real x1,x2,x3,x4
      integer jt(pcols)         
      integer maxg(pcols)      
      integer ideep(pcols)   
      integer lengath       
c
      real, parameter :: pmin = 1.e-20   ! minimum precip rate (ms-1)
      real, parameter :: epsq = 2.e-12
      real, parameter :: rcp  = 1.0/cpair
      real, parameter :: cmf0 = 100./gravit
      real, parameter :: a0   = 0.200 !coefficient of eq(4) from slingos(1991) (default 0.245)
      real, parameter :: b0   = 0.120 !coefficient of eq(4) from slingos(1991) (default 0.125)
      real, parameter :: cvmax = 0.75 !upper limit of convectice cloud cover   (default 0.8)
      real, parameter :: elocp = latvap/cpair  ! for evap. of cloud water 
c
c---------------------------local workspace-----------------------------
c
c most variables with appended 2 are a second copy of similar quantities
c declared above to accommodate calls to two convection schemes
c
      integer  i,k                ! lon, lev, constituent indices
c
c-----------------------------------------------------------------------
      dt2 = 2.0 * delt
      rdt = 1.0/dt2
c
c reset
c
      detr(1:imx22,1:kmx)=0.
      clw(1:imx22,1:kmx)=0.
      rn(1:imx22)=0.
c
c zero out precip and convective fields before accumulating terms
c
      precc(1:pcols) = 0.
      snow(1:pcols) = 0.
      ptend%s(1:pcols,1:pver)=0.
      ptend%q(1:pcols,1:pver,1)=0.
      cmfdqr(1:pcols,1:pver)=0.
      cmfmc(1:pcols,1:pver)=0.
c
c!!! note the ccmcnv is top down, whereas ncep is bottom up 
c
      psfc(1:imx2) = ps(1:imx2) * 1000.  ! [cb]->[Pa]

      state%pint(1:imx2,1) = 1.0         ! avoid log(pint(i,1))=log(0)
      do k = 2, pver
         kk =  pver - k + 2
         do i = 1, imx2
           state%pint(i,k) = psfc(i) * si(i,kk)
         enddo
      enddo
      state%pint(1:imx2,pverp) = psfc(1:imx2)

      do k = 1, pver
         do i = 1, imx2
           state%pmid(i,k) = 0.5*(state%pint(i,k)+state%pint(i,k+1))
         enddo
      enddo
      do k = 1, pver
         kk = pver - k + 1
         do i = 1, imx2
           state%t(i,k)   = t1(i,kk)
           state%q(i,k,1) = q1(i,kk)
           ptendt(i,k)    = (t1(i,kk)-tp(i,kk))*rdt
           ptendq(i,k)    = (q1(i,kk)-qp(i,kk))*rdt
         enddo
      enddo
c
c set layer thicknesses
c
      do k = 1,pver
         do i = 1,imx2
            state%pdel(i,k) = state%pint(i,k+1) - state%pint(i,k)
         end do
      end do
      pmln(1:imx2,1:pver)  = log(state%pmid(1:imx2,1:pver))
      piln(1:imx2,1:pverp) = log(state%pint(1:imx2,1:pverp))
c
      state%phis(1:imx2)  = 0.      ! assume zero terrain
      state%lchnk = lat
      state%ncol  = imx2
c
c calculate geopotential height 
c
      call virtem(state%t ,state%q ,zvir    ,tv      )
      call zmid(piln(1,pverp) ,pmln    ,rair    ,gravit   ,tv  ,
     $          state%zm      )
      call zint(piln    ,pmln    ,rair    ,gravit  ,tv      ,
     $          state%zm      ,state%zi      ) 
c
c Initialization of ZM constants
c
      call zm_convi ( x1, x2, x3, x4 )
c
c Main of ZM
c
      call zm_convr(state%lchnk   ,state%ncol   ,
     1  state%t ,state%q    ,precc      ,cnt         ,cnb     ,
     2  hpbl    ,state%zm   ,state%phis ,state%zi    ,ptend%q ,
     3  ptend%s ,state%pmid ,state%pint ,state%pdel  ,
     4  delt    ,cmfmc   ,cmfdqr  ,
     5  gamt    ,dlf     ,pflx    ,zdu     ,rprd    ,
     6  mu      ,md      ,du      ,eu      ,ed      ,
     7  dp      ,dsubcld ,jt      ,maxg    ,ideep   ,
     8  lengath ,ql      ,rliq    ,ptendq  ,ptendt )

      cmfmc(1:imx2,1:pver)=cmfmc(1:imx2,1:pver)*cmf0  ! [mb/s]->[kg/m^2/s]

#ifdef EVAP
c
c  update
c
      state%t(1:imx2,1:pver) = state%t(1:imx2,1:pver)
     1           + ptend%s(1:imx2,1:pver) * rcp * dt2
      state%q(1:imx2,1:pver,1) = state%q(1:imx2,1:pver,1)
     1           + ptend%q(1:imx2,1:pver,1)       * dt2

      ptend%s(1:pcols,1:pver)=0.   ! reset
      ptend%q(1:pcols,1:pver,1)=0. ! reset
c
c     cldfrc is calculated like cnvcldn.F based on slingos(1991)
c
      do i = 1, imx2
        if (precc(i).gt.0.0) then
          xx=precc(i)*86400.e3         ! [m/s] -> [mm/day]
          xx=a0+b0*alog(xx)            ! convective cloud fraction 
          cldfrc(i,1:pver)=max(0.,min(xx,cvmax)) ! same for all levels
        else
          cldfrc(i,1:pver)=0.0
        endif
      end do
!
! Evaporate some of the precip directly into the environment (Sundqvist)
      call zm_conv_evap(state, ptend, rprd, cldfrc, dt2,
     1                  precc, snow, .false.)
#endif

c
c output
c
      kbot(1:imx2) = pver - nint(cnb(1:imx2)) + 1
      ktop(1:imx2) = pver - nint(cnt(1:imx2)) + 1
      do i = 1,imx2
         if(precc(i).gt.pmin.and.kbot(i).lt.ktop(i)) then
            rn(i) = precc(i) * dt2   ! [m/s] -> [m]
            kuo(i) = 1
         else
            rn(i) = 0.
            kuo(i) = 0
         endif
      end do
      do k = 1, pver
         kk = pver - k + 1
         do i = 1, imx2
           if(kuo(i).eq.1) then
             t1(i,k) = state%t(i,kk) + ptend%s(i,kk) * rcp * dt2
             q1(i,k) = state%q(i,kk,1) + ptend%q(i,kk,1)       * dt2
             detr(i,k) = zdu(i,kk)*dt2  ! [-]
             clw(i,k)  = dlf(i,kk)*dt2  ! [kg/kg]
           endif
         enddo
      enddo

#ifdef CLD1
c
c treatment of detrained cloud water
c   -> evaporation for all.
c      if oversaturated, it will fall as grid-scale rain.
c
      do k = 1, pver
        do i = 1, imx2
          if (clw(i,k).gt.0.0) then
            q1(i,k) = q1(i,k) + clw(i,k)
            t1(i,k) = t1(i,k) - elocp*clw(i,k)
            clw(i,k) = 0.0
          endif
        end do
      end do
#endif
c
c----------------store t, q  for next time step
c
      do k = 1, kmx
        do i = 1, imx2
          tp(i,k)  = tp1(i,k)
          qp(i,k)  = qp1(i,k)
          tp1(i,k) = t1(i,k)
          qp1(i,k) = max(q1(i,k),epsq)
        enddo
      enddo

#ifdef DBG
      i = max(imx2/4,1)
      write(6,*)
      write(6,*) '      precip/hr ',           precc(i) * 3600.
      write(6,*)
      write(6,*)' kbot ktop kuo ', kbot(i),ktop(i),kuo(i)
 144  format(a,2g10.3,a)
      write(6,*)
      write(6,*)
      write(6,145)
 145  format('lev' ,8x,'pres',4x,' t  ',7x,'  q  ',5x,'mass flx',
     $     3x,' t  tend',2x,' q  tend')
      write(6,146)
 146  format('   ' ,8x,' pa ',4x,' k  ',7x,'kg/kg',5x,' kg/m^2s',
     $     3x,' deg/s  ',2x,' kg/kg/s')
      write(6,*)
      do k=1,levs_
         write(6,155)k,state%pmid(i,k),state%t(i,k),
     $        state%q(i,k,1),cmfmc(i,k),
     $        ptend%s(i,k)*rcp,ptend%q(i,k,1)
      enddo
 155  format(i2,5x,f9.2,1x,f8.2,5(2x,g10.3),3x,a)
#endif
      return
      end
 
