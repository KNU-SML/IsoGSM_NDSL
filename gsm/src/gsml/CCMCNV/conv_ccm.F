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
      subroutine conv_ccm(t       ,qh      ,pcpc    ,jctop   ,jcbot   ,
     $                    pblh    ,zm      ,geos    ,zi      ,qtg     ,
     $                    ttg     ,pap     ,paph    ,dpp     ,ts      ,
     $                    delt    ,mcon    ,cme     ,nstep   ,lat     ,
     $                    tpert   ,qpert   ,dlf     ,pflx    ,cape)
#include <paramodel.h>
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
C
C ##### MAIN DRIVER FOR ZHANG-MCFARLANE CONVECTION SCHEME #####
C
C-----------------------------------------------------------------------
C This is contributed code not fully standardized by the CCM core group.
C All variables have been typed, where most are identified in comments
C The current procedure will be reimplemented in a subsequent version 
C of the CCM where it will include a more straightforward formulation 
C and will make use of the standard CCM nomenclature
C-----------------------------------------------------------------------
C
C $Id: conv_ccm.F,v 1.2 2009/10/02 00:35:46 kana Exp $
C
#include <implicit.h>
C******************************************************************

C same as conv.up except saturation vapor pressure is calculated
C in a different way.
C
C jul 17/92 - guang jun zhang, m.lazare. calls new buoyan, q1q2
C             and moment (several new work fields added for later).
C
C nov 21/91 - m.lazare. like previous conv except calls new
C                       clpdprp.
C feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.
C             previous version conv.
C performs deep convective adjustment based on mass-flux closure
C algorithm.
C
C ************************ index of variables **********************
C  i      => input arrays.
C  i/o    => input/output arrays.
C  w      => work arrays.
C  wg     => work arrays operating only on gathered points.
C  ic     => input data constants.
C  c      => data constants pertaining to subroutine itself.
C
C  wg * alpha    array of vertical differencing used (=1. for upstream).
C  wg * betad    downward mass flux at cloud base.
C  wg * betau    upward   mass flux at cloud base.
C  w  * cape     convective available potential energy.
C  wg * capeg    gathered convective available potential energy.
C  c  * capelmt  threshold value for cape for deep convection.
C  ic  * cpres    specific heat at constant pressure in j/kg-degk.
C  i  * dpp      local sigma half-level thickness (i.e. dshj).
C  ic  * delt     length of model time-step in seconds.
C  wg * dp       layer thickness in mbs (between upper/lower interface).
C  wg * dqdt     mixing ratio tendency at gathered points.
C  wg * dsdt     dry static energy ("temp") tendency at gathered points.
C  wg * dudt     u-wind tendency at gathered points.
C  wg * dvdt     v-wind tendency at gathered points.
C  wg * dsubcld  layer thickness in mbs between lcl and maxi.
C  ic  * grav     acceleration due to gravity in m/sec2.
C  wg * du       detrainment in updraft. specified in mid-layer
C  wg * ed       entrainment in downdraft.
C  wg * eu       entrainment in updraft.
C  wg * hmn      moist static energy.
C  wg * hsat     saturated moist static energy.
C  w  * ideep    holds position of gathered points vs longitude index.
C  ic  * levs_     number of model levels.
C  ic  * ilg      lon+2 = size of grid slice.
C  wg * j0       detrainment initiation level index.
C  wg * jd       downdraft   initiation level index.
C  ic  * jlat     gaussian latitude index.
C  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
C  wg * jt       top  level index of deep cumulus convection.
C  ic  * kount    current model timestep number.
C  w  * lcl      base level index of deep cumulus convection.
C  wg * lclg     gathered values of lcl.
C  w  * lel      index of highest theoretical convective plume.
C  wg * lelg     gathered values of lel.
C  w  * lon      index of onset level for deep convection.
C  wg * long     gathered values of lon.
C  ic  * lev      levs_+1.
C  w  * maxi     index of level with largest moist static energy.
C  wg * maxg     gathered values of maxi.
C  wg * mb       cloud base mass flux.
C  wg * mc       net upward (scaled by mb) cloud mass flux.
C  wg * md       downward cloud mass flux (positive up).
C  wg * mu       upward   cloud mass flux (positive up). specified 
C                at interface
C  ic  * msg      number of missing moisture levels at the top of model.
C  c  * kups     number of points undergoing deep convection.
C  w  * p        grid slice of ambient mid-layer pressure in mbs.
C  i  * pblt     row of pbl top indices.
C  i/o * pcp      row of precipitable water in metres.
C  w  * pcpdh    scaled surface pressure.
C  w  * pf       grid slice of ambient interface pressure in mbs.
C  wg * pg       grid slice of gathered values of p.
C  i  * pressg   row of surface pressure in pa.
C  w  * q        grid slice of mixing ratio.
C  wg * qd       grid slice of mixing ratio in downdraft.
C  wg * qdb      row of qd at cloud base.
C  wg * qg       grid slice of gathered values of q.
C  i/o * qh       grid slice of specific humidity.
C  w  * qh0      grid slice of initial specific humidity.
C  wg * qhat     grid slice of upper interface mixing ratio.
C  wg * ql       grid slice of cloud liquid water.
C  wg * qs       grid slice of saturation mixing ratio.
C  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
C  wg * qstpg    grid slice of gathered values of qstp.
C  wg * qu       grid slice of mixing ratio in updraft.
C  ic  * rgas     dry air gas constant.
C  wg * rl       latent heat of vaporization.
C  w  * s        grid slice of scaled dry static energy (t+gz/cp).
C  wg * sd       grid slice of dry static energy in downdraft.
C  wg * sdb      row of sd at cloud base.
C  wg * sg       grid slice of gathered values of s.
C  wg * shat     grid slice of upper interface dry static energy.
C  i  * shbj     grid slice of local bottom interface sigma values.
C  i  * shj      grid slice of local half-level sigma values.
C  i  * shtj     row of local top interfaces of first level.
C  wg * su       grid slice of dry static energy in updraft.
C  wg * sumde    row of vertically-integrated moist static energy 
C                change.
C  wg * sumdq    row of vertically-integrated scaled mixing ratio 
C                change.
C  wg * sumdt    row of vertically-integrated dry static energy change.
C  wg * sumq     row of vertically-integrated mixing ratio change.
C  i/o * t        grid slice of temperature at mid-layer.
C  o  * jctop    row of top-of-deep-convection indices passed out.
C  o  * jcbot    row of base of cloud indices passed out.
C  w  * tf       grid slice of temperature at interface.
C  wg * tg       grid slice of gathered values of t.
C  w  * tl       row of parcel temperature at lcl.
C  wg * tlg      grid slice of gathered values of tl.
C  w  * tp       grid slice of parcel temperatures.
C  wg * tpg      grid slice of gathered values of tp.
C  i/o * u        grid slice of u-wind (real).
C  wg * ug       grid slice of gathered values of u.
C  i/o * utg      grid slice of u-wind tendency (real).
C  i/o * v        grid slice of v-wind (real).
C  w  * va       work array re-used by called subroutines.
C  wg * vg       grid slice of gathered values of v.
C  i/o * vtg      grid slice of v-wind tendency (real).
C  i  * w        grid slice of diagnosed large-scale vertical velocity.
C  w  * z        grid slice of ambient mid-layer height in metres.
C  w  * zf       grid slice of ambient interface height in metres.
C  wg * zfg      grid slice of gathered values of zf.
C  wg * zg       grid slice of gathered values of z.
C
C-----------------------------------------------------------------------
C
C multi-level i/o fields:
C
C input/output arguments:
C
      real t(ILOTS,levs_) 
      real qh(ILOTS,levs_) 
c      real utg(ILOTS,levs_) 
c      real vtg(ILOTS,levs_) 
      real qtg(ILOTS,levs_) 
      real ttg(ILOTS,levs_)
C
C       input arguments
      real pap(ILOTS,levs_) 
      real paph(ILOTS,levs_+1) 
      real dpp(ILOTS,levs_) 
      real zm(ILOTS,levs_) 
      real geos(ILOTS) 
      real zi(ILOTS,levs_+1)
      real pblh(ILOTS) 
      real zs(ILOTS) 
      real tpert(ILOTS) 
      real qpert(ILOTS)
C
C       output arguments
      real pcpck(ILOTS,levs_)
      real mup(ILOTS,levs_) 
      real mdn(ILOTS,levs_) 
      real mcon(ILOTS,levs_) 
      real dlg(ILOTS,levs_)    ! gathered version of the detraining cld h2o tend
      real dsdt2(ILOTS,levs_)  ! scattered version of the temp tend
      real dqdt2(ILOTS,levs_)  ! scattered version of the q tend
      real dlg2(ILOTS,levs_)   ! gathered version of the detraining cld h2o tend
      real dlf(ILOTS,levs_)    ! scattered version of the detraining cld h2o tend
      real pflx(ILOTS,levs_)   ! scattered precip flux at each level
      real pflxg(ILOTS,levs_)  ! gather precip flux at each level
      real cu(ILOTS,levs_)     ! scattered condensation rate
      real cug(ILOTS,levs_)    ! gathered condensation rate 
      real evpg(ILOTS,levs_)   ! gathered evap rate of rain in downdraft
      real evp(ILOTS,levs_)    ! scattered evap rate of rain in downdraft
      real mumax(ILOTS) 
      real cme(ILOTS,levs_)
      real small, msn, ms(ILOTS)
      integer km1, kp1
C
C single-level i/o fields:
C
C       input arguments
      real ts(ILOTS) 
      real pblt(ILOTS)
C
C       input/output arguments:
C
      real paprc(ILOTS) 
      real paprs(ILOTS) 
      real zctopm(ILOTS)
C
C       output arguments:
C
      real jctop(ILOTS) 
      real jcbot(ILOTS)
      real pcpr(ILOTS) 
      real pcps(ILOTS) 
      real pcpc(ILOTS)
C
C-----------------------------------------------------------------------
C
C general work fields (local variables):
C
      real u(ILOTS,levs_) 
      real v(ILOTS,levs_) 
      real q(ILOTS,levs_) 
      real p(ILOTS,levs_) 
      real z(ILOTS,levs_) 
      real s(ILOTS,levs_) 
      real qh0(ILOTS,levs_) 
      real tp(ILOTS,levs_) 
      real zf(ILOTS,levs_+1) 
      real pf(ILOTS,levs_+1) 
      real qstp(ILOTS,levs_) 

      real cape(ILOTS) 
      real tl(ILOTS) 
      real sumq(ILOTS) 
      real pcpdh(ILOTS)
c      real sumdt(ILOTS) 
c      real sumdq(ILOTS) 
c      real sumde(ILOTS)

      integer lcl(ILOTS) 
      integer lel(ILOTS) 
      integer lon(ILOTS) 
      integer maxi(ILOTS) 
      integer ideep(ILOTS) 
      integer index(ILOTS)
      real precip
C
C gathered work fields:
C
      real qg(ILOTS,levs_) 
      real tg(ILOTS,levs_) 
      real pg(ILOTS,levs_) 
      real zg(ILOTS,levs_) 
      real sg(ILOTS,levs_) 
      real tpg(ILOTS,levs_) 
      real zfg(ILOTS,levs_+1) 
      real qstpg(ILOTS,levs_) 
      real ug(ILOTS,levs_) 
      real vg(ILOTS,levs_) 
      real cmeg(ILOTS,levs_)

      real capeg(ILOTS) 
      real tlg(ILOTS)

      integer lclg(ILOTS) 
      integer lelg(ILOTS) 
      integer maxg(ILOTS)
C
C work fields arising from gathered calculations.
C
      real mu(ILOTS,levs_) 
      real eu(ILOTS,levs_) 
      real dqdt(ILOTS,levs_) 
      real dsdt(ILOTS,levs_) 
      real du(ILOTS,levs_) 
      real md(ILOTS,levs_) 
      real ed(ILOTS,levs_) 
      real alpha(ILOTS,levs_) 
      real sd(ILOTS,levs_) 
      real qd(ILOTS,levs_) 
      real mc(ILOTS,levs_) 
      real qhat(ILOTS,levs_) 
      real qu(ILOTS,levs_) 
      real su(ILOTS,levs_) 
      real qs(ILOTS,levs_) 
      real shat(ILOTS,levs_) 
      real dp(ILOTS,levs_) 
      real hmn(ILOTS,levs_) 
      real hsat(ILOTS,levs_) 
      real ql(ILOTS,levs_) 
      real dudt(ILOTS,levs_) 
      real dvdt(ILOTS,levs_) 
      real ud(ILOTS,levs_) 
      real vd(ILOTS,levs_)

c      real deltat(ILOTS,levs_) 
c      real deltaq(ILOTS,levs_)

      real betau(ILOTS) 
      real betad(ILOTS) 
      real qdb(ILOTS) 
      real sdb(ILOTS) 
      real dsubcld(ILOTS) 
      real mb(ILOTS) 
      real totpcp(ILOTS) 
      real totevp(ILOTS)
      real mu2(ILOTS,levs_) 
      real eu2(ILOTS,levs_) 
      real du2(ILOTS,levs_) 
      real md2(ILOTS,levs_) 
      real ed2(ILOTS,levs_)
      real denom
      integer jt(ILOTS) 
      integer jlcl(ILOTS) 
      integer j0(ILOTS) 
      integer jd(ILOTS)
      real capelmt
      real cpres
      real delt
      real fixer
      integer i
      integer ii
      integer k
      integer lat
      integer lengath
      integer msg
      integer nstep
      real psdiss
      real psevap
      real psheat
      real psrain
      real qdifr
      real qeff
      real qmin
      real rl
      real sdifr
C
C--------------------------Commons--------------------------------------
C 
#include <guang.h>
#include <comcmf.h>
C
C--------------------------Data statements------------------------------
C
      logical momentm
      data capelmt/70./
      data momentm/.FALSE./
C
C set up convection cap, this was moved from mfinti to here to
C get rid of dependency on hybrid pressure coefficients for the
C Suarez repos code
C
      if (pap(1,1) .ge. 4.e3) then 
         limcnv =1
      else
        do k=1,levs_-1
          if (pap(1,k).lt.4.e3 .and. pap(1,k+1).ge.4.e3) then
            limcnv = k
            goto 1010
          end if
        end do
        limcnv = levs_+1
      end if
 1010 continue
c
c Set internal variable "msg" (convection limit) to "limcnv-1"
c
      msg = limcnv - 1
C
C initialize necessary arrays.
C zero out variables not used in ccm
C
      do i = 1,ILOTS
         paprc(i) = 0.
         paprs(i) = 0.
         zctopm(i) = 0.
      end do
      psdiss = 0.
      psheat = 0.
      psevap = 0.
      psrain = 0.
c      jlatpr = 32
      cpres = 1004.64
      a = 21.656
      b = 5418.
      c1 = 6.112
      c2 = 17.67
      c3 = 243.5
      eps1 = 0.622
      qmin = 1.E-20
      tfreez = 273.16
      rl = 2.5104E6
C
C initialize convective tendencies
C
      do k = 1,levs_
         do i = 1,ILOTS
            dqdt(i,k) = 0.
            dsdt(i,k) = 0.
            dudt(i,k) = 0.
            dvdt(i,k) = 0.
c            deltaq(i,k) = qh(i,k)
c            deltat(i,k) = t(i,k)
            dqdt2(i,k) = 0.
            dsdt2(i,k) = 0.
            pcpck(i,k) = 0.
            pflx(i,k) = 0.
            pflxg(i,k) = 0.
         end do
      end do
      if (.not.momentm) then
         do k = msg + 1,levs_
            do i = 1,ILOTS
               u(i,k) = 0.
               v(i,k) = 0.
            end do
         end do
      end if
C
      do i = 1,ILOTS
         pblt(i) = levs_
         pcpr(i) = 0.
         pcps(i) = 0.
         dsubcld(i) = 0.
         sumq(i) = 0.
c         sumdt(i) = 0.
c         sumdq(i) = 0.
         pcpdh(i) = rgrav
         jctop(i) = levs_
         jcbot(i) = 1
      end do
C
C calculate local pressure (mbs) and height (m) for both interface
C and mid-layer locations.
C
      do i = 1,ILOTS
         zs(i) = geos(i)*rgrav
         pf(i,levs_+1) = paph(i,levs_+1)*0.01
         zf(i,levs_+1) = zi(i,levs_+1) + zs(i)
      end do
      do k = 1,levs_
         do i = 1,ILOTS
            p(i,k) = pap(i,k)*0.01
            pf(i,k) = paph(i,k)*0.01
            z(i,k) = zm(i,k) + zs(i)
            zf(i,k) = zi(i,k) + zs(i)
         end do
      end do
C
      do k = levs_ - 1,msg + 1,-1
         do i = 1,ILOTS
            if (abs(z(i,k)-zs(i)-pblh(i)).lt.
     $          (zf(i,k)-zf(i,k+1))*0.5) pblt(i) = k
         end do
      end do
C
C store incoming specific humidity field for subsequent calculation
C of precipitation (through change in storage).
C convert from specific humidity (bounded by qmin) to mixing ratio.
C define dry static energy (normalized by cp).
C
      do k = 1,levs_
         do i = 1,ILOTS
            qh0(i,k) = qh(i,k)
            qeff = max(qh(i,k),qmin)
            q(i,k) = qeff
            s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
            tp(i,k)=0.0
            shat(i,k) = s(i,k)
            qhat(i,k) = q(i,k)
            dp(i,k) = dpp(i,k)*0.01
            qg(i,k) = q(i,k)
            tg(i,k) = t(i,k)
            pg(i,k) = p(i,k)
            zg(i,k) = z(i,k)
            sg(i,k) = s(i,k)
            tpg(i,k) = tp(i,k)
            zfg(i,k) = zf(i,k)
            qstpg(i,k) = q(i,k)
            ug(i,k) = u(i,k)
            vg(i,k) = v(i,k)
            dlg(i,k) = 0
            dlg2(i,k) = 0
            dlf(i,k) = 0
            cu(i,k) = 0
            evp(i,k) = 0
         end do
      end do
      do i = 1,ILOTS
         zfg(i,levs_+1) = zf(i,levs_+1)
         capeg(i) = 0.
         lclg(i) = 1
         lelg(i) = levs_
         maxg(i) = 1
         tlg(i) = 400.
         dsubcld(i) = 0.
         qdb(i) = 0.
         sdb(i) = 0.
         betau(i) = 0.
         betad(i) = 0.
      end do
C
C evaluate covective available potential energy (cape).
C
      call buoyan(q       ,t       ,p       ,z       ,pf       ,
     $            tp      ,qstp    ,tl      ,rl      ,cape     ,
     $            pblt    ,lcl     ,lel     ,lon     ,maxi     ,
     $            rgas    ,grav    ,cpres   ,msg     ,nstep    ,
     $            lat     ,tpert   ,qpert)
C
C determine whether grid points will undergo some deep convection
C (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
C (require cape.gt. 0 and lel<lcl as minimum conditions).
C
         call whenfgt(ILOTS,cape,1,capelmt,index,lengath)
         if (lengath.eq.0) return
CDIR$ IVDEP
         do ii=1,lengath
            i=index(ii)
            ideep(ii)=i
         end do
cc        jyes = 0
cc        jno = ILOTS - 1 + 2
cc        do il = 1,ILOTS
cc           if (cape(il).gt.capelmt) then
cc              jyes = jyes + 1
cc              ideep(jyes) = il
cc           else
cc              jno = jno - 1
cc              ideep(jno) = il
cc           end if
cc        end do
cc        lengath = jyes
cc        if (lengath.eq.0) return
C
C obtain gathered arrays necessary for ensuing calculations.
C
      do k = 1,levs_
         do i = 1,lengath
            dp(i,k) = 0.01*dpp(ideep(i),k)
            qg(i,k) = q(ideep(i),k)
            tg(i,k) = t(ideep(i),k)
            pg(i,k) = p(ideep(i),k)
            zg(i,k) = z(ideep(i),k)
            sg(i,k) = s(ideep(i),k)
            tpg(i,k) = tp(ideep(i),k)
            zfg(i,k) = zf(ideep(i),k)
            qstpg(i,k) = qstp(ideep(i),k)
            ug(i,k) = u(ideep(i),k)
            vg(i,k) = v(ideep(i),k)
         end do
      end do
C
      do i = 1,lengath
         zfg(i,levs_+1) = zf(ideep(i),levs_+1)
      end do
      do i = 1,lengath
         capeg(i) = cape(ideep(i))
         lclg(i) = lcl(ideep(i))
         lelg(i) = lel(ideep(i))
         maxg(i) = maxi(ideep(i))
         tlg(i) = tl(ideep(i))
      end do
C
C calculate sub-cloud layer pressure "thickness" for use in
C closure and tendency routines.
C
      do k = msg + 1,levs_
         do i = 1,lengath
            if (k.ge.maxg(i)) then
               dsubcld(i) = dsubcld(i) + dp(i,k)
            end if
         end do
      end do
C
C define array of factors (alpha) which defines interfacial
C values, as well as interfacial values for (q,s) used in
C subsequent routines.
C
      do k = msg + 2,levs_
         do i = 1,lengath
            alpha(i,k) = 0.5
            sdifr = 0.
            qdifr = 0.
            if (sg(i,k).gt.0. .or. sg(i,k-1).gt. 0.) 
     $          sdifr = abs((sg(i,k)-sg(i,k-1))/
     $          max(sg(i,k-1),sg(i,k)))
            if (qg(i,k).gt.0. .or. qg(i,k-1).gt.0.) 
     $          qdifr = abs((qg(i,k)-qg(i,k-1))/
     $          max(qg(i,k-1),qg(i,k)))
            if (sdifr.gt.1.E-6) then
               shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/
     $                     (sg(i,k-1)-sg(i,k))
            else
               shat(i,k) = 0.5* (sg(i,k)+sg(i,k-1))
            end if
            if (qdifr.gt.1.E-6) then
               qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/
     $                     (qg(i,k-1)-qg(i,k))
            else
               qhat(i,k) = 0.5* (qg(i,k)+qg(i,k-1))
            end if
         end do
      end do
C
C obtain cloud properties.
C
      call cldprpccm(qg      ,tg      ,ug      ,vg      ,pg      ,
     $            zg      ,sg      ,mu      ,eu      ,du      ,
     $            md      ,ed      ,sd      ,qd      ,ud      ,
     $            vd      ,mc      ,qu      ,su      ,zfg     ,
     $            qs      ,hmn     ,hsat    ,alpha   ,shat    ,
     $            ql      ,totpcp  ,totevp  ,cmeg    ,maxg    ,
     $            lelg    ,jt      ,jlcl    ,maxg    ,j0      ,
     $            jd      ,rl      ,1       ,lengath ,rgas    ,
     $            grav    ,cpres   ,msg     ,nstep   ,lat     ,
     $            pflxg   ,evpg    ,cug     ,mu2     ,eu2     ,
     $            du2     ,md2     ,ed2     ) 
C
C determine cloud base mass flux.
C
      do i = 1,lengath
         qdb(i) = qd(i,maxg(i))
         sdb(i) = sd(i,maxg(i))
         betad(i) = md(i,maxg(i))
         betau(i) = mu(i,maxg(i))
      end do
C
C convert detrainment from units of "1/m" to "1/mb".
C
      do k = msg + 1,levs_
         do i = 1,lengath
            du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
            ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         end do
      end do

      call closure(qg      ,tg      ,pg      ,zg      ,sg      ,
     $             tpg     ,qs      ,qu      ,su      ,mc      ,
     $             du      ,mu      ,md      ,qd      ,sd      ,
     $             alpha   ,qhat    ,shat    ,dp      ,qstpg   ,
     $             zfg     ,ql      ,dsubcld ,mb      ,capeg   ,
     $             tlg     ,lclg    ,lelg    ,jt      ,maxg    ,
     $             1       ,lengath ,rgas    ,grav    ,cpres   ,
     $             rl      ,msg     ,capelmt ,nstep   ,lat     )
C
C limit cloud base mass flux to theoretical upper bound.
C
      if (.true.) then
         do k = msg + 2,levs_
            do i = 1,lengath
               if (mu(i,k).gt.0.) then
                  mb(i) = min(mb(i),dp(i,k)/ (2.*delt*mu(i,k)))
               end if
            end do
         end do
      else

         do i = 1,lengath
            mumax(i) = 0
         end do

         do k = msg + 2,levs_
            do i = 1,lengath
               mumax(i)  = max(mumax(i), mu(i,k)/dp(i,k))
            end do
         end do

         do i = 1,lengath
            if (mumax(i).gt.0.) then
               mb(i) = min(mb(i),0.5/(delt*mumax(i)))
            else
               mb(i) = 0.
            endif
         end do
      endif

      do k = msg + 1,levs_
         do i = 1,lengath
            mu(i,k) = mu(i,k)*mb(i)
            md(i,k) = md(i,k)*mb(i)
            mc(i,k) = mc(i,k)*mb(i)
            du(i,k) = du(i,k)*mb(i)
            eu(i,k) = eu(i,k)*mb(i)
            ed(i,k) = ed(i,k)*mb(i)
            cmeg(i,k) = cmeg(i,k)*mb(i)
            cug(i,k) = cug(i,k)*mb(i)
            evpg(i,k) = evpg(i,k)*mb(i)
            pflxg(i,k) = pflxg(i,k)*mb(i)*100./grav
            mu2(i,k) = mu2(i,k)*mb(i)
            md2(i,k) = md2(i,k)*mb(i)
            du2(i,k) = du2(i,k)*mb(i)
            eu2(i,k) = eu2(i,k)*mb(i)
            ed2(i,k) = ed2(i,k)*mb(i)
         end do
      end do
#ifdef DBG
c perform some consistency checks on updraft quantities
      small = 1.e-18
      do i = 1,lengath
         ms(i) = 0
      end do
      do k = levs_,1,-1
         kp1 = min(k+1,levs_)
         do i = 1,lengath
            msn = ms(i) + (eu2(i,k)-du2(i,k))*dp(i,k)
            if (abs(msn-mu2(i,k))/(abs(msn)+abs(mu2(i,k))+small)
     $           .gt.1.e-10.and.mu2(i,k).ne.0.) then
               write (6,*) ' updraft consistency ', i, k, lat, msn, 
     $              ms(i), mu2(i,k), mu2(i,kp1), eu2(i,k), du2(i,k),
     $              (zfg(i,k)-zfg(i,k+1)), dp(i,k)
               write (6,*) ' dmu ', msn-ms(i), mu2(i,k)-mu2(i,kp1)
               write (6,*) ' jb or maxg ', maxg(i)
            endif
            ms(i) = msn
         end do
      end do
c perform some consistency checks on downdraft quantities
      do i = 1,lengath
         ms(i) = 0
      end do
      do k = 1,levs_
         kp1 = min(k+1,levs_)
         do i = 1,lengath
            msn = ms(i) - (ed2(i,k))*dp(i,k)
            denom = max(abs(msn),abs(md2(i,kp1)),abs(ed2(i,k)*dp(i,k)))
            if (     abs(msn-md2(i,kp1)).gt.1.e-10*denom
     $          .and.denom.gt.1.e-15
     $         ) then
               write (6,*) ' inconsistent downdraft ', i, k, lat, msn, 
     $              ms(i), md2(i,kp1), md2(i,k), ed2(i,k)
               write (6,*) ' jt, jd, jb ', jt(i), jd(i), maxg(i)
               write (6,*) ' rdiff is ',
     $          abs(msn-md2(i,kp1))/(abs(msn)+abs(md2(i,kp1))+1.e-16)
            endif
            ms(i) = msn
         end do
      end do
c      if (lat.eq.32) then
c         stop
c      endif
C
#endif
      do i = 1,lengath
         betau(i) = betau(i)*mb(i)
         betad(i) = betad(i)*mb(i)
C
C totpcp from cldprp has the dimension of kg/kg, here it is 
C converted to kg/(m^2*s), the precipitation rate
C
         totpcp(i) = totpcp(i)*mb(i)*100./grav
         totevp(i) = totevp(i)*mb(i)*100./grav
      end do
C
C compute temperature and moisture changes due to convection.
C
      call q1q2(dqdt    ,dsdt    ,qg      ,sg      ,qs      ,
     $          qu      ,su      ,mc      ,du      ,alpha   ,
     $          qhat    ,shat    ,dp      ,mu      ,md      ,
     $          sd      ,qd      ,ql      ,dsubcld ,qdb     ,
     $          sdb     ,betau   ,betad   ,mb      ,lclg    ,
     $          jt      ,maxg    ,2.*delt ,1       ,lengath ,
     $          cpres   ,rl      ,msg     ,nstep   ,lat     ,
     $          dlg     ,cug     )
C
C compute momentum changes due to convection, if desired (i.e
C if logical switch set).
C
c
c      if(momentm)                                                   
c       then
c        call moment(dudt,dvdt,du,alpha,dp,ed,eu,mc,md,mu,
c     1             pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,
c     2             dsubcld,maxg,jd,jt,rl,
c     3             msg,2.*delt,grav,cpres,rgas,levs_,1,lengath,ILOTS,lat)
c      endif
c
c      if (.false.) then
      if (ntrac_.gt.1) then
         call convtran
     $               (qh      ,mu2     ,md2     ,
     $                du2     ,eu2     ,ed2     ,dp      ,
     $                dsubcld ,jt      ,maxg    ,ideep   ,
     $                1       ,lengath ,nstep   ,lat     ,
     $                delt    )
      endif
C
C gather back temperature and mixing ratio.
C
      do k = msg + 1,levs_
         do i = 1,lengath
            psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+
     $               dvdt(i,k)*v(ideep(i),k))*dpp(ideep(i),k)/grav
C
C q is updated to compute net precip, and then reset to old value.
C the last line is overwritten. so the input basic variables, i.e.
C q, t, u and v are updated to include convective increments. 
C (5/2/95)
C
            q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
            t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)
            u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
            v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
            cme(ideep(i),k) = cmeg(i,k)
            mup(ideep(i),k) = mu(i,k)
            mdn(ideep(i),k) = md(i,k)
            mcon(ideep(i),k) = mc(i,k)
            qtg(ideep(i),k) = dqdt(i,k)
            ttg(ideep(i),k) = dsdt(i,k)
c            utg(ideep(i),k) = dudt(i,k)
c            vtg(ideep(i),k) = dvdt(i,k)
            dlf(ideep(i),k) = dlg(i,k)
            pflx(ideep(i),k) = pflxg(i,k)
            cu(ideep(i),k) = cug(i,k)
            evp(ideep(i),k) = evpg(i,k)
         end do
      end do
C
      do i = 1,lengath
         jctop(ideep(i)) = jt(i)
         jcbot(ideep(i)) = jlcl(i)
         psevap = psevap + totevp(i)
         psrain = psrain + totpcp(i)
      end do
C
C convert back to specific humidity from mixing ratio.
C take into account any moisture added to ensure positiveness
C of specific humidity at start of routine.
C
      do k = msg + 1,levs_
         do i = 1,ILOTS
            qh(i,k) = q(i,k)
            qh(i,k) = qh(i,k) - max((qmin-qh0(i,k)),0.)
         end do
      end do
C
C determine change in column storage due to deep convection.
C enforce water vapor conservation if negative precip diagnosed
C
      do k = levs_,msg + 1,-1
         do i = 1,ILOTS
            sumq(i) = sumq(i) - dpp(i,k)* (qh(i,k)-qh0(i,k))
            if (sumq(i).lt.0.0) then
               fixer = -sumq(i)/dpp(i,k)
               qh(i,k) = qh(i,k) - fixer
               t(i,k) = t(i,k) + (rl/cpres)*fixer
               sumq(i) = 0.0
C
               qtg(i,k) = qtg(i,k) - fixer/ (2.0*delt)
               ttg(i,k) = ttg(i,k) + (rl/cpres)*fixer/ (2.0*delt)
            end if
            pcpck(i,k) = max(0.,sumq(i))
         end do
      end do
C
C obtain final precipitation rate.
C
      do i = 1,ILOTS
c         llo1 = ts(i) .ge. tfreez
C
C here pcpr and pcps are in units of kg/m^2, ie. precip per
C time step
C
c         pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),0.),0.,llo1)
c         pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),0.),llo1)
         precip = pcpdh(i)*max(sumq(i),0.)
         if (ts(i) .ge. tfreez) then
           pcpr(i) = precip
           pcps(i) = 0.
         else
           pcpr(i) = 0.
           pcps(i) = precip
         end if
      end do

C
C accumulate precipitation, the 1000. is the density of water, so
C paprc and paprs are now in units of meters.
C
      do i = 1,ILOTS
         paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
         paprs(i) = paprs(i) + pcps(i)/1000.
      end do
C
C convert precipitation to m/s, ie, precip rate.
C
      do i = 1,ILOTS
         pcpr(i) = pcpr(i)/ (2.*delt)/1000.
         pcps(i) = pcps(i)/ (2.*delt)/1000.
         pcpc(i) = pcpr(i) + pcps(i)
         psheat = psheat + (pcps(i)+pcpr(i))*rl
      end do
      do k = msg + 1,levs_
         do i = 1,ILOTS
            pcpck(i,k) = pcpdh(i)*pcpck(i,k)/ (2.*delt)
         end do
      end do
C
C calculate conservation of quantities.
C

C       if(lat.eq.jlatpr)then
c        do l=msg+1,levs_
c        do i=1,lengath
c          sumdq(i) = sumdq(i) + 2.*delt*(rl/cpres)*dpp(ideep(i),l)*
c     1                          dqdt(i,l)
c          sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
c        end do
c        end do
c
c        write(6,*)'sumdq,sumdt,sumde in convection subroutine########'
c        do i=1,lengath
c          sumde(i) = sumdt(i) + sumdq(i)
c          write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
c        end do
cc
c        write(6,*)'sumdq,sumdt,sumde ... all points'
c      do i=1,ILOTS
c         sumdq(i) = 0.0
c         sumdt(i) = 0.0
c      end do
cc
c      do l=msg+1,levs_
c      do i=1,ILOTS
c        deltaq(i,l) = qh(i,l) - deltaq(i,l)
c        deltat(i,l) = t (i,l) - deltat(i,l)
c        sumdq(i) = sumdq(i) + (rl/cpres)*dpp(i,l)*deltaq(i,l)
c        sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
c      end do
c      end do
c      do i=1,ILOTS
c        sumde(i) = sumdt(i) + sumdq(i)
c      end do
c        write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,ILOTS)
c
c  901   format(1x,3e20.12, i10, i10)
c  902   format(1x,i10, 3e20.12)
c
c      endif

      return
      end
 
