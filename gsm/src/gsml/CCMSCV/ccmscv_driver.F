#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NTRAC ntrac
      subroutine ccmscv_driver(imx2,imx22,kmx,
     &              jcap,delt,si,del,sl,slk,ps,nstep,
     &              q1,t1,rn,kbot,ktop,kuo,spd,lat,slimsk,dot,
     &              hpbl,gamt,gamq,qci,qrs,delx)
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#define ILOTS igrd12p_
#else
#define ILOTS lonf2p_
#endif
#else
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
c-----------------------------------------------------------------------
c   ncar ccm3 mid/shallow convection
c-----------------------------------------------------------------------
#include <implicit.h>
c------------------------------parameters-------------------------------
#include <comadj.h>
c------------------------------arguments--------------------------------
c
c input arguments
c
      integer nstep                 ! current time step for cmfmca diag.
      integer lat                   ! latitude index (s->n)
      real tv(ILOTS,levs_)           ! virtual temperature array
      real zm(ILOTS,levs_)           ! height above sfc in mid-layer
      real phis(ILOTS)               ! surface geopotential
      real pmid(ILOTS,levs_)         ! pressures at model levels
      real pint(ILOTS,levs_+1)       ! pressure at model interfaces
      real pdel(ILOTS,levs_)         ! delta-p (pressure depth)
      real rpdel(ILOTS,levs_)        ! 1./pdel
      real pmln(ILOTS,levs_)         ! log(pmid)
      real piln(ILOTS,levs_+1)        ! log(pint)
      real psfc(ILOTS)                ! surface temperature
c
c input/output arguments
c
      real t(ILOTS,levs_)            ! temperature (k)
      real q(ILOTS,levs_,ntrac_)     ! constituent mixing ratio field
      real ts(ILOTS)                ! surface temperature
c
c output arguments
c
      real precc(ILOTS)             ! convective-scale preciptn rate
      real cnt(ILOTS)               ! top level of convective activity
      real cnb(ILOTS)               ! bottom level of convective activity
      real cmfdt(ILOTS,levs_)        ! dt/dt due to moist convection
      real cmfdq(ILOTS,levs_)        ! dq/dt due to moist convection
      real zmdt(ILOTS,levs_)         ! zhang convective temperature tendency
      real zmdq(ILOTS,levs_)         ! zhang convective moisture tendency
      real cmfdqr(ILOTS,levs_)       ! dq/dt due to moist convective rainout 
      real cmfmc(ILOTS,levs_)        ! moist convection cloud mass flux
      real cmfsl(ILOTS,levs_)        ! moist convection lw stat energy flux
      real qc(ILOTS,levs_)            ! dq/dt due to rainout terms
      real cmflq(ILOTS,levs_)         ! Moist convection total water flux
      real tdt
c
      real del(imx22,kmx),sl(imx22,kmx),slk(imx22,kmx),
     $          ps(imx2),si(imx22,kmx+1),
     &          q1(imx22,kmx),t1(imx22,kmx),rn(imx2)
      integer   kbot(imx2),ktop(imx2),kuo(imx2)
      real spd(imx2,kmx),slimsk(imx2),dot(imx2,kmx)
      real hpbl(imx2),gamt(imx2),gamq(imx2)
      real qci(imx2,kmx),qrs(imx2,kmx)
c
      integer   kbots(ILOTS),ktops(ILOTS),kuos(ILOTS)
      data pmin /1.e-20/   ! minimum precip rate (ms-1)
c
c---------------------------local workspace-----------------------------
c
c most variables with appended 2 are a second copy of similar quantities
c declared above to accommodate calls to two convection schemes
c
      integer  i,k                ! lon, lev, constituent indices
c
c-----------------------------------------------------------------------
c
c define thermodynamic constansts
c
      call initcommons
c
c zero out precip and convective fields before accumulating terms
c
      call resetr(precc,ILOTS,0.)
      call resetr(cmfdt ,ILOTS*levs_,0.)
      call resetr(cmfdq ,ILOTS*levs_,0.)
      call resetr(zmdt  ,ILOTS*levs_,0.)
      call resetr(zmdq  ,ILOTS*levs_,0.)
      call resetr(cmfdqr,ILOTS*levs_,0.)
      call resetr(cmfmc ,ILOTS*levs_,0.)
      call resetr(cmfsl ,ILOTS*levs_,0.)
      call resetr(qc ,ILOTS*levs_,0.)
      call resetr(cmflq ,ILOTS*levs_,0.)
c
      tdt = delt * 2.0
c
c!!! note the ccmcnv is top down, whereas ncep is bottom up 
c
c
c   convert surface pressure to pascal from cb
c
      do i = 1, ILOTS
         psfc(i) = ps(i) * 1000.
         pint(i,levs_+1) = psfc(i)
      enddo
      do k = 1, levs_
         kk =  levs_ - k + 2
         do i = 1, ILOTS
           pint(i,k) = psfc(i) * si(i,kk)
         enddo
      enddo
      do k = 1, levs_
         do i = 1, ILOTS
           pmid(i,k) = 0.5*(pint(i,k)+pint(i,k+1))
         enddo
      enddo
      do k = 1, levs_
         kk = levs_ - k + 1
         do i = 1, ILOTS
           t(i,k) = t1(i,kk)
           q(i,k,1) = q1(i,kk)
         enddo
      enddo
c
c set layer thicknesses
c
      do k = 1,levs_
         do i = 1,ILOTS
            pdel(i,k) = pint(i,k+1) - pint(i,k)
         end do
      end do
      do k = 1,levs_
         do i = 1,ILOTS
            pmln(i,k) = log(pmid(i,k))
            rpdel(i,k) = 1./pdel(i,k)
            piln(i,k) = log(pint(i,k))
         end do
      end do
      do i = 1,ILOTS
         piln(i,levs_+1) = log(pint(i,levs_+1))
      end do
c
      do i = 1,ILOTS
         phis(i) = 0.  ! assume zero terrain
         ts(i) = 400.  ! no check for snow or rain
      end do
c
c calculate geopotential height 
c
      call virtem(t       ,q       ,zvir    ,tv      )
      call zmid(piln(1,levs_+1),pmln    ,rair    ,gravit   ,tv  ,
     $          zm      )
      call cmfmca(lat     ,nstep   ,tdt    ,pmid    ,pdel    ,
     $            rpdel   ,zm      ,gamt   ,gamq   ,phis    ,
     $            hpbl   ,t       ,q       ,cmfdt   ,cmfdq   ,
     $            cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   ,
     $            qc      ,cnt     ,cnb    )
c
c output
c
      do i = 1,ILOTS
         kbots(i) = levs_ - nint(cnb(i)) + 1
         ktops(i) = levs_ - nint(cnt(i)) + 1
      end do
      do i = 1,ILOTS
         if(kbots(i).lt.ktops(i)) then
            rn(i) = precc(i) * tdt         ! m
            kuos(i) = 1
         else
            rn(i) = 0.
            kuos(i) = 0
         endif
      end do
c
c  add precipitating convection to radiation
c
      do i = 1,ILOTS
         if(kuos(i).eq.1.and.rn(i).gt.pmin) then
            kuo(i) = 1
            kbot(i) = min(kbot(i),kbots(i))
            ktop(i) = max(ktop(i),ktops(i))
         endif
      end do
c
c  largescale feedback of nonprecip(precip) clouds 
c
      do k = 1, levs_
         kk = levs_ - k + 1
         do i = 1, ILOTS
           if(kuos(i).eq.1) then
              t1(i,k) = t(i,kk)
              q1(i,k) = q(i,kk,1)
           endif
         enddo
      enddo
#ifdef DBG
      i = ILOTS/4
      write(6,*)
      write(6,*) ' precip/h ', precc(i)*3600
      write(6,*)
      write(6,*)' kbot ktop kuo ', kbots(i),ktops(i),kuos(i)
 144  format(a,2g10.3,a)
      write(6,*)
      write(6,*)
      write(6,145)
 145  format('lev' ,8x,'pres',4x,' t  ',7x,'  q  ',5x,'mass flx',
     $     3x,' t  tend',2x,' q  tend')
      write(6,146)
 146  format('   ' ,8x,' pa ',4x,' k  ',7x,'kg/kg',5x,' kg/m^2s',
     $     3x,' deg/s  ',2x,' kg/m^2s')
      write(6,*)
      do k=1,levs_
         write(6,155)k,pmid(i,k),t(i,k),q(i,k,1),cmfmc(i,k),
     $        cmfdt(i,k),cmfdq(i,k)
      enddo
 155  format(i2,5x,f9.2,1x,f8.2,5(2x,g10.3),3x,a)
#endif
      return
      end
 
