#define G g
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define RD rd
      subroutine gwdrag(ims2,imx2,kmx,
     .                  u1,v1,t1,pstar,vtj,usqj,ksm,kbj,kbps,
     .                  velco,bnvl2,roll,ro,taub,si,del,sl,sigk,rcl,
     .                  lat,kdt,hprimx,akwnmb,tensio)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gwdrag      performs gravity wave drag computations.
c   prgmmr: jordan c. alpert org: w/nmc23    date: 91-03-12
c
c abstract: performs gravity wave drag computations.
c
c program history log:
c   87-06-03  jordan c. alpert  gwdrag, gwdps.
c   89-02-01  hann-ming henry juang  change to fortran 77.
c   91-03-03  sela-rozwodoski  cray guard code and constants -
c                              ggwdra, ggwdps.
c   02-09-09  yifeng for rsm mp
c
c usage:    call gwdrag(u1,v1,t1,pstar,vtj,usqj,ksm,kbj,kbps,
c                       velco,bnvl2,roll,ro,taub,si,del,sl,sigk,rcl,
c                       lat,kdt,hprimx,akwnmb,tensio)
c   input argument list:
c     u1       - zonal      wind component *cos(lat)  m/sec at t0-dt.
c     v1       - meridional wind component *cos(lat)  m/sec at t0-dt.
c     t1       - temperature deg k at t0-dt.
c     pstar    - surface pressure (cb).
c     vtj      - virtual temperature.
c     usqj     - richardson number.
c     ksm      - top    of low level layer for gwdrag set to 2.
c     kbj      - bottom of low level layer for gwdrag set to 1.
c     kbps     - bottom starting sigma level for p&s stress calculation
c     si(n)    - p/psfc at base of layer n.
c     del(n)   - positive increment of p/psfc across layer n.
c     sl(n)    - p/psfc at middle of layer n.
c     rcl      - reciprocal of square of cos(lat).
c     lat      - latitude  number.  used only as diagnostic.
c     kdt      - time step number.  used only as diagnostic.
c     hprimx   - topographic standard deviation  (m).
c
c   output argument list:
c     velco    - comp of wind along the direction of low level layer.
c     bnvl2    - brunt-viasila freq as (n) and also (n**2).
c     roll     - low level ro.
c     ro       - density  mts.
c     taub     - surface stress.
c     sigk     - single dimension array of length kdim which holds
c              - the constants of inverse sigma values raised to
c              - r/cp power from subroutine ggwdps.
c     akwnmb   - length scale.
c     tensio   - stress.
c
c remarks: list caveats, other helpful hints or information
c
c attributes:
c   language: fortran, cft77.
c   machine:  cray y-mp.
c
c$$$
c
c  g l a s   m i g w d  (gwdrag)
c
#include <paramodel.h>
      dimension
     .   u1(imx2,kmx),v1(imx2,kmx),t1(imx2,kmx),
     .   pstar(imx2),vtj(imx2,kmx),usqj(imx2,kmx),
     .   velco(imx2,kmx-1),bnvl2(imx2),roll(imx2),
     .   ro(imx2,kmx),taub(imx2),
     .   si(imx2,kmx+1),del(imx2,kmx),sl(imx2,kmx),sigk(imx2,kmx),
     .   hprimx(imx2),tensio(imx2,kmx)
c local array
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
      dimension hprime(ILOTS)
     .,         hco(ILOTS),hsi(ILOTS)
     .,         crif2(ILOTS),fro2(ILOTS)
     .,         cl(levs_),
     .          bnv2(ILOTS),ulow(ILOTS)
     .,         zmean(levs_,2)
c
      logical icrilv(ILOTS)
      real akappa
c hmhj
      idimt=ims2
      idimt2=imx2
      kdim=kmx
c
c-----only do constants first time through monnin
c
c
c-----gsfc to nmc bridge constants
c
      ntopm1=kdim-1
c
        grav = g_
        grav2 = grav * grav
        agrav = 1./grav
        rgas = rd_
        gr2 = 2. * grav2 / rgas
        akappa = 2. / 7.
        frocut=0.85 * 0.85
c
        rcs=sqrt(rcl)
c
c-----the variance of topography (comes in as std dev) on a lat pair
c
      do 100 i= 1, idimt
      hprime(i) = hprimx(i) * hprimx(i)
 100  continue
c
c-----constrain variance to be not greater than 160000 m**2
c     only if p&s low level stress is used - comment for gfcl low lev
c
      do 200 i=1, idimt
      value = 1.6 e 5
      hprime(i) = min( hprime(i), value )
 200  continue
c
c-----initialize critical level control vector bits all to zero
c
       do 300 i=1, idimt
          icrilv(i) = .false.
 300   continue
c
c----bnvl2, the low level brunt-viasla frequency is not a fct of k
c----sqrt (n**2) nb: i use bnvl2 for n**2 and now n itself
c
        do 400 i=1,idimt
        value = 0.0
        bnvl2(i) = max( bnvl2(i) , value )
        bnvl2(i) = sqrt ( bnvl2(i) )
 400    continue
c
c-----set initial values for stress
c
        do 500 k=1,kdim
        do 500 i=1,idimt
          tensio(i,k) = 0.0
 500    continue
c
c------ level loop
c
c-----set up bottom values of stress if we are not starting from
c-----from level 1 - nb (-) here by convention
c
       do 600 klow = 1, kbps
       do 600 i=1,idimt
          tensio(i,klow) = -1. * taub(i)
 600   continue
c
          do 1100 k = kbps, ntopm1
c
          do 1110 i=1, idimt
          fro2(i) = 0.0
cra 1110     continue
c
c-----calculate squared brunt vaisala frequency at level k
c---- sigkm is recip( sigma** kappa) at top of layer
c---- sigkin is recip( sigma** kappa) at bottom of layer
c
          sigkm = sigk(i,k+1)
          sigkin = sigk(i,k)
c
c--- n**2 as function of k - branch on low level <ksm - tv is used (vtj)
c or just allow this type of operation on usqj  (ri)
c
cra         do 1120 i=1,idimt
            bnv2(i) = gr2 * (sl(i,k+1) + sl(i,k))
     .    * (vtj(i,k+1) * sigkm - vtj(i,k) * sigkin)
     .    /  (  (vtj(i,k+1) * sigkm
     .                 + vtj(i,k) * sigkin)
     .       * (sl(i,k) - sl(i,k+1))
     .                 * (t1(i,k+1) + t1(i,k)) )
c
c-----unstable layer if ri < 0 - - using ri in place of n**2
c
        icrilv(i) = icrilv(i) .or. ( usqj(i,k)
     1                    .le. 0.0 )
c
       icrilv(i) = icrilv(i) .or. ( usqj(i,k)
     1                                  .lt. 0.25 )
cra 1120   continue
c
c-----compute critical froude    ---> means stable (icrilv=0)
c
cra         do 1150 i=1, idimt
          if    ( .not. icrilv(i) )  then
                crif2(i) = 1 - .25 / usqj(i,k)
                crif2(i) = crif2(i) * crif2(i)
          else
          crif2(i)=0.0
          endif
cra 1150    continue
c
c-----unstable layer if upper air vel comp along surf vel <=0 (crit lay)
c---- at (u-c)=0. crit layer exists and bit vector should be set (.le.)
c
cra        do 1160 i=1,idimt
       icrilv(i) = icrilv(i) .or. velco(i,k)
     1                                  .le. 0.0
cra 1160  continue
c
c-----sqrt (n**2)  nb: i use bnv2 for n**2 and now n  itself
c
cra        do 1170 i=1, idimt
         value = 0.
         bnv2(i) = max( bnv2(i) , value )
         bnv2(i) = sqrt ( bnv2(i) )
cra 1170   continue
 1110   continue
c
c-----computing stress at surface  and 1 level up & limit max value
c
c----using taub at kbj level
c
       if(k .eq. kbj)  then
c
        do 1180 i=1, idimt
         crif2(i) = min( crif2(i) , frocut )
cra 1180   continue
cra        do 1190 i=1, idimt
               if    ( .not. icrilv(i) )  then
c
c-----gfdl low level surface stress is negitive with respect to glas
c-----becasue glas subtracts tendency while monin adds  we change taub
c-----to -taub (nb but change back for monin tendency
c
                 tensio(i,k) = -1. * taub(i)
c
c---fr**2 at surface only
c
               fro2(i) = bnvl2(i) * bnvl2(i) *
     1        hprime(i) / (velco(i,k) * velco(i,k))
c
               endif
cra 1190   continue
 1180   continue
c
       else
c
c-----in glas version there is a calculation at the boundary layer
c-----which is not a "part" of the model std levels and a calc at the
c-----first model layer.  the value for tensio at the bndy layer is set
c-----to the value at the first model layer.  in nmc the bndy layer can
c-----be ignored all together because the stress (tensio) is on sis,
c-----interfaces while the deacceration is on layers sls.
c
c-----compute the local froude   for the stable case - make sure
c-----that by chance the projection of the local wind (u1,v1)
c-----should not be smaller then 1.m/s, since it is cubed in the
c-----denominator
c-----
c
          do 1195 i=1,idimt
            if   ( .not.  icrilv(i) ) then
                 fro2(i) = bnv2(i) / ( (akwnmb * 0.5) *
     1          ( ro(i,k) + ro(i,k+1) ) *
     2      velco(i,k) * velco(i,k) * velco(i,k) ) *
     3          tensio(i,k)
            endif
 1195     continue
c
       endif
c
c-----compute stress at level in question for stable case
c
         do 1210 i=1,idimt
            if( .not.icrilv(i) .and. fro2(i) .gt. crif2(i) ) then
c
c------fro2 changed to> from .ge.
c
                  tensio(i,k+1) = tensio(i,k) * crif2(i) / fro2(i)
            endif
c
c----- constant stress if crit froude   not met (.le. 6/1)
c
            if( .not. icrilv(i) .and. fro2(i) .le. crif2(i) )  then
                  tensio(i,k+1) = tensio(i,k)
             endif
 1210    continue
c
c-----all done - pass back stress profile and vertically diff
c
 1100 continue
c
      return
      end
