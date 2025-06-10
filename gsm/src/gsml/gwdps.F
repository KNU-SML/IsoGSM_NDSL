#define CP cp
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
#define RV rv
      subroutine gwdps(ims2,imx2,kmx,a,b,
     1                  u1,v1,t1,q1,
     2                  pstar,
     3                  si,del,cl,sl,rcl,deltim,lat,kdt,hprime,
     4                  dusfc,dvsfc)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gwdps       includes gravity wave drag.
c   prgmmr: jordan c. alpert org: w/nmc23    date: 91-03-12
c
c abstract: using the gwd parameterizations of ps-glas and ph-
c   gfdl technique, the time tendencies of u v
c   are altered to include the effect of mountain induced
c   gravity wave drag from sub-grid scale orography including
c   convective breaking, shear breaking and the presence of
c   critical levels.
c
c program history log:
c   87-06-03  jordan c. alpert  fr30(v3h-mx)
c   02-09-09  yifeng cui for rsm mp
c
c usage:    call gwdps(a,b,u1,v1,t1,q1,pstar,
c                      si,del,cl,sl,rcl,deltim,lat,kdt,hprime)
c   input argument list:
c     a        - negative non-lin tendency for v wind component.
c     b        -          non-lin tendency for u wind component.
c     u1       - zonal      wind component *cos(lat)  m/sec at t0-dt.
c     v1       - meridional wind component *cos(lat)  m/sec at t0-dt.
c     t1       - temperature deg k at t0-dt.
c     q1       - specific humidity at t0-dt.
c     pstar    - surface pressure (cb).
c     si(n)    - p/psfc at base of layer n.
c     del(n)   - positive increment of p/psfc across layer n.
c     cl(n)    = 1 - sl(n).
c     sl(n)    - p/psfc at middle of layer n.
c     rcl      - reciprocal of square of cos(lat).
c     deltim   - time step  secs.
c     lat      - latitude  number.
c     kdt      - time step number.
c     hprime   - topographic standard deviation  (m).
c
c   output argument list:
c     a        - as augmented by tendency due to migwd.
c     b        - as augmented by tendency due to migwd.
c
c   output files:
c     ft06f001 - printout file.
c
c attributes:
c   language: fortran, cft77.
c   machine:  cray y-mp.
c
c$$$
c
#include <paramodel.h>
      dimension a(imx2,kmx),  b(imx2,kmx),
     2         u1(imx2,kmx), v1(imx2,kmx),hprime(imx2),
     3         t1(imx2,kmx), q1(imx2,kmx),pstar(imx2),
     7          si(imx2,kmx+1),del(imx2,kmx),cl(imx2,kmx),sl(imx2,kmx)
      dimension dusfc(imx2),dvsfc(imx2)
      dimension delks(imx2),delks1(imx2)
c
      parameter(cp=cp_,g=g_,rd=rd_,rv=rv_)
      parameter(gor=g/rd,gocp=g/cp,fv=rv/rd-1)
      parameter(dw2min=1.,rimin=-100.)
c
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
      logical ldrag(ILOTS)
      dimension rdzt(ILOTS,levs_-1)
      dimension velco(ILOTS,levs_-1)
c
c  debug dimension
c
      dimension ri(levs_-1)
c
c----         mountain induced gravity wave drag
c----   unit14  - subgrid scale mountain variance height input
c----   common to be added to smf,gloo for migwd        *j*
c----   idimt = 256 and kdim = 18 for example...
c
      dimension taub(ILOTS),xn(ILOTS),
     *          yn(ILOTS),
     *          dtaux(ILOTS,levs_),
     *          ubar(ILOTS),vbar(ILOTS),
     *          taud(ILOTS,levs_),
     *          fr(ILOTS),gf(ILOTS),
     *          ulow(ILOTS),bnv(ILOTS),
     *          vtj(ILOTS,levs_),
     *          dtauy(ILOTS,levs_),
     *          bnv2(ILOTS,levs_),
     *          sigk(ILOTS,levs_),
     *          taup(ILOTS,levs_+1),
     *          usqj(ILOTS,levs_),
     *          ro(ILOTS,levs_),
     *          roll(ILOTS)
      dimension dtfac(ILOTS)
      dimension velko(levs_-1)
      real akappa
      real slreal
c hmhj
      idimt=ims2
      idimt2=imx2
      kdim=kmx
c
c
c--------------------------------------------constants for migwd *j*
c
c-----only do constants first time through monnin
c
      rlowlv =0.7
      xl     =4.0 e 4
      critac =5.0 e -4
      ncnt   =100
c
            akappa = 2. / 7.
c
      do j=1,idimt
      do 18 k=1,kdim
      slreal = sl(j,k)
      sigk(j,k)  = slreal**(-akappa)
   18 continue
      enddo
c
c------- mkdimp   the top sigma level over which migwd will operate
c
       mkdimp = kdim+1
c
c------- kbj is the bottom of the low 1/3 level usually = 1
c
       ksm=1
       do j=1,idimt2
        kbj = 1
         do 15 k = kbj, kdim
c         if ((si(j,k) .lt. rlowlv).and.(ksm.lt.k)) then
c         ksm = k
c         go to 16
c         endif
         if ((si(j,k) .ge. rlowlv).and.(k.gt.ksm-1)) ksm=k+1
   15    continue
c   16   continue
       enddo
c       print*,'ksm',ksm
c
c-----ksm -1 intervals in the lower third of atm (sigma < .667)
c
       ksmm1 = ksm - 1
       do j=1,idimt
         delks(j) = si(j,kbj)-si(j,ksm)
         delks1(j)= sl(j,kbj)-sl(j,ksm)
       enddo
c
c----above, the low layer delta sigma
c-----below the starting sigma level for ps stress calc defaults to 2
c
       kbps = 2
       lcap=kdim
       lcapp1 = lcap + 1
       factop=0.5
c
        grav = g
        grav2 = grav * grav
        rgas = rd
        gr2 = 2.0 * grav2 / rgas
        gmax = 1.
        aj   = 1.
        xlinv =  1.0 / xl
        veleps=1.0
      rcs=sqrt(rcl)
      cs = 1. / rcs
c
c----------saving richardson number in usqj for migwd        *j*
c
      do k=1,kdim-1
        do j=1,idimt
          rdzt(j,k)=gor*si(j,k+1)/(sl(j,k)-sl(j,k+1))
        enddo
      enddo
      do k=1,kdim
        do j=1,idimt
          vtj(j,k)=t1(j,k)*(1.+fv*q1(j,k))
        enddo
      enddo
      do k=1,kdim-1
        do j=1,idimt
          ti=0.5*(t1(j,k)+t1(j,k+1))
          rdz=rdzt(j,k)/ti
          dw2=rcl*((u1(j,k)-u1(j,k+1))**2+(v1(j,k)-v1(j,k+1))**2)
          shr2=max(dw2,dw2min)*rdz**2
          bvf2=g*(gocp+rdz*(vtj(j,k+1)-vtj(j,k)))/ti
          usqj(j,k)=max(bvf2/shr2,rimin)
        enddo
      enddo
c
c-----vertical structure of ri in ri(kdim) for diagnostics
c
c-----the linear mountain induced gravity mode p&s prameterization
c      explitly done
c-----this routine computes the deceleration of the zonal wind and
c-----meridional wind due to mountain gravity drag.
c
c-----    code variables          description
c
c-----       xn,yn            projections of "low-level" wind
c-----                        in zonal & meridional directions
c
c-----       ulow             "low-level" wind magnitude -        (= u)
c-----                        averaged up to 2km above surface
c
c-----       bnv2             bnv2 = n**2
c
c-----       hprime           sub-grid scale mountain height      (= h)
c-----                        from navy tape, averaged,'envelope'std. va
c-----                        read in in smf,common-ed to gloo
c
c-----       taub             base momentum flux
c-----                        = -(ro * u**3/(n*xl)*gf(fr) for n**2> 0
c-----                        = 0.                        for n**2 < 0
c
c
c-----       fr               froude    =   n*hprime / u
c-----       g                gmax*fr**2/(fr**2+aj**2)
c-----       gmax             = 1.0
c-----       aj               = 1.0
c
c-----ksm is defined as the number of levels up 1/3 from the lowest used
c-----to calculate the "low-level" averages.
c
c
c-----initialize arrays     (on cyber)
c
       do 200 i=1,idimt
         xn(i) = 0.0
         yn(i) = 0.0
       ubar (i)=  0.0
       vbar (i)=  0.0
       roll (i)=  0.0
       taub (i)=  0.0
       ulow (i)=  0.0
       taup (i,kdim+1) = 0.0
 200   continue
c
       do 250 k=1,kdim
       do 250 i=1,idimt
        taup(i,k) = 0.0
        ro(i,k) = sl(i,k) * pstar(i)
     1                 / ( rgas * vtj(i,k) )
 250  continue
c
c----density   tons/meter**3
c--------.---------.---------.---------.---------.---------.---------. .
c-----compute low level averages
c-----(u,v)*cos(lat)  use uv=(u1,v1) which is wind at t0-1
c----- use rcs=1/cos(lat) to get wind field
c----  ksm   the top of the lowest 1/3 layer "the low level" is 6
c
        do 300 k=kbj,ksmm1
        do 300 i=1,idimt
        rcsks = rcs * del(i,k) / delks(i)
        ubar(i) = ubar(i) + rcsks * u1(i,k)
        vbar(i) = vbar(i) + rcsks * v1(i,k)
 300    continue
c
c----compute the "low level" or 1/3 wind magnitude (m/s)
c
      do 400 i=1,idimt
        ulow(i) = sqrt( ubar(i) * ubar(i) + vbar(i) * vbar(i) )
 400  continue
      do 450 i=1,idimt
        value = 1.0
        ulow(i) = max( ulow(i), value )
 450  continue
c
c-----calculate squared low level brunt vaisala frequency over the
c-----first ksm levels then average
c---- sigkm is recip( sigma** kappa) at top of lower layer
c---- sigkin is recip( sigma** kappa) at bottom of lower layer
c---- rdelks (del(k)/delks) vert ave factor so we can * instead of /
c
        do 500 i=1,idimt
           bnv2(i,1) = 0.
 500    continue
c
        do 550 k=kbj,ksmm1
        do 550 i=1,idimt
            bnv2(i,k) = gr2 * (sl(i,k) + sl(i,k+1))
     1    * (vtj(i,k+1) * sigk(i,k+1) - vtj(i,k) * sigk(i,k))
     2    /  (  (vtj(i,k) * sigk(i,k)
     3                 + vtj(i,k+1) * sigk(i,k+1))
     4       * (sl(i,k) - sl(i,k+1))
     5                 * (t1(i,k) + t1(i,k+1)) )
 550    continue
c
        do 600 k=1,kdim-1
        do 600 i=1,idimt
        velco(i,k) =
     1    (0.5*rcs)*( (u1(i,k) + u1(i,k+1)) *
     2                      ubar(i) +
     3                       (v1(i,k) + v1(i,k+1)) *
     4                      vbar(i))
        velco(i,k)=velco(i,k)/ulow(i)
        if ((velco(i,k).lt.veleps).and.(velco(i,k).ge.0.)) then
          velco(i,k) = veleps
        endif
 600    continue
c
c  no drag when critical level in the base layer
c
        do 700 i=1,idimt
        ldrag(i)=velco(i,1).le.0.
 700    continue
        do 750 k=2,ksmm1
        do 750 i=1,idimt
        ldrag(i)=ldrag(i).or. velco(i,k).le.0.
 750    continue
c
c  no drag when bnv2.lt.0
c
        do 800 k=1,ksmm1
        do 800 i=1,idimt
        ldrag(i)=ldrag(i).or. bnv2(i,k).lt.0.
 800    continue
c
c-----the low level weighted average ri is stored in usqj(1,1; idimt)
c-----the low level weighted average n**2 is stored in bnv2(1,1; idimt)
c----- this is called bnvl2 in gwdrag not bnv2
c
        kbjp1 = kbj + 1
        do 900 i=1,idimt
        wtkbj = (sl(i,kbj)-sl(i,kbjp1))/delks1(i)
        usqj(i,1) = wtkbj * usqj(i,kbj)
        bnv2(i,1) = wtkbj * bnv2(i,kbj)
 900    continue
c
        do 1000 k = kbjp1,ksmm1
        do 1000 i=1,idimt
           rdelks = (sl(i,k)-sl(i,k+1))/delks1(i)
           bnv2(i,1) =   bnv2(i,1) + bnv2(i,k)  * rdelks
           usqj(i,1) =   usqj(i,1) + usqj(i,k)   * rdelks
 1000   continue
c
        do 1010 i=1,idimt
        ldrag(i)=ldrag(i).or. bnv2(i,1).le.0.0
        ldrag(i)=ldrag(i).or. ulow(i).eq.1.0
 1010   continue
c
c ----- set all ri low level values to the low level value
c
        kbjbeg = kbj
        if(kbj .eq. 1) kbjbeg = 2
        do 1020 k=kbjbeg,ksmm1
        do 1020 i=1,idimt
           usqj(i,k) = usqj(i,1)
 1020   continue
c
c-----  low level density
c
        do 1030 k=kbj,ksmm1
        do 1030 i=1,idimt
           rdelks =  del(i,k)/delks(i)
           roll(i) =   roll(i) + ro(i,k)   * rdelks
 1030   continue
c
        do 1050 i=1,idimt
         if  (.not.ldrag(i) )   then
c
c-----vector square root function - vsqrt -  used to compute bnv
c
       bnv(i) = sqrt( bnv2(i,1) )
c
c-----calculate fr  froude    ---- n*hprime / u
c
       fr(i) = bnv(i) * hprime(i) / ulow(i)
c
c----continue w/ where block
c
c-----calculate g   the universal flux function
c
       gf(i) = gmax * fr(i) * fr(i) /
     1                 ( fr(i) * fr(i) + aj * aj )
c
c-----calculate taub - (the base flux)
c-----remember - the low level n is in bnv2(1,1;idimt) = bnv = bnvl2
c
       taub(i) =  -xlinv * roll(i) *
     1     ulow(i) * ulow(i) * ulow(i) * gf(i) / bnv(i)
c
c-----calculate xn, yn
c
               xn(i) = ubar(i) / ulow(i)
c
               yn(i) = vbar(i) / ulow(i)
c
          else
c
              taub(i) = 0.0
               xn(i)  = 0.0
               yn(i)  = 0.0
c
          endif
 1050   continue
c
c------the call to gwdrag:
c------      taup are returned    other parameters from monn
c
      call gwdrag(idimt,idimt2,kdim,
     1            u1,v1,t1,pstar,vtj,usqj,ksm,kbj,kbps,
     1            velco,bnv2,roll,ro,taub,si,del,sl,sigk,rcl,
     2            lat,kdt,hprime,xlinv,taup)
c
           if(lcap.lt.kdim) then
c
           do 1100 klcap =  lcapp1, kdim
           do 1100 i=1,idimt
           sira =          si(i,klcap) / si(i,lcap)
           taup(i,klcap) =  sira * taup(i,lcap)
c
 1100      continue
c
           endif
c
c-----fix up the level 1 (or more) stress to be linear with level
c-----kbps and the stress at the bottom taub (if kbps is .gt. 1)
c
        if (kbps .gt. kbj) then
           kbpsm1 = kbps - 1
           kbpsp1 = kbps + 1
c
           do 1200 ik1 = kbj, kbpsm1
           do 1200 i=1,idimt
           savem   =  ( (si(i,ik1+1) - si(i,kbpsp1)) /
     1                (si(i,kbpsp1) - 1.0 ) )
           taup(i,ik1+1) =  taup(i,kbpsp1)   -
     1        savem * ( taup(i,kbj) - taup(i,kbpsp1) )
 1200      continue
c
        endif
c
c---keep in mind that taup is zero-ed out before each call
c-----vertically difference stress for d tau / d sigma  from si to sl
c
      do 1300 k=1,kdim
      do 1300 i=1,idimt
c
c-----the stress in gwsdrag has been calc using -taub which now must be
c-----returned to -(amplitude) below the old way for tau as -(ro*u**3/nl
c-----instead of ro*uamp*k*n*hprime**2
c
c
      taud(i,k) = (taup(i,k+1) - taup(i,k) ) / del(i,k)
c
c---where del= si(k)-si(k+1)  (sign 'switched' in sela code)
c
 1300 continue
c
c-----calculate deceleration terms - dtaux,dtauy
c
      do 1400 k=1,kdim
      do 1400 i=1,idimt
      taud(i,k) = taud(i,k) / pstar(i)
 1400 continue
c
c------limit de-acceleration (momentum deposition ) at top to 1/2 value
c------the idea is some stuff must go out the 'top'
c
c------limit de-acceleration (momentum deposition ) at top to 1/2 value
c------the idea is some stuff must go out the 'top'
c
      do 1500 klcap = lcap, kdim
      do 1500 i=1,idimt
      taud(i,klcap) = taud(i,klcap) * factop
 1500 continue
c
c----- *g and * by cos(lat) for mrf tendencies
c
      csgrav = cs * grav
      do 1600 k=1,kdim
      do 1600 i=1,idimt
      taud(i,k) = taud(i,k) * csgrav
 1600 continue
c
c------if the gravity wave drag would force a critical line
c------in the lower ksmm1 layers during the next 2*deltim timestep,
c------then only apply drag until that critical line is reached.
c
      do 1610 i=1,idimt
      dtfac(i)=1.
1610  continue
      do 1620 k=1,ksmm1
      do 1620 i=1,idimt
      if(taud(i,k).ne.0.)
     &dtfac(i)=min(dtfac(i),abs(velco(i,k)/(2.*deltim*rcs*taud(i,k))))
1620  continue
1625  continue
      do 1630 k=1,kdim
      do 1630 i=1,idimt
      taud(i,k)=taud(i,k)*dtfac(i)
1630  continue
c
           do 1660 k=1,kdim
           do 1660 i=1,idimt
c
            dtaux(i,k) = xn(i) * taud(i,k)
            dtauy(i,k) = yn(i) * taud(i,k)
c
 1660      continue
c
c-----done with calculation - add it to old a and old b
c-----a corresponds to dtauy term and b to dtaux
c
      do 2001 i=1,idimt
      dusfc(i)=0.
      dvsfc(i)=0.
 2001 continue
      do 2000 k=1,kdim
      do 2000 i=1,idimt
      a(i,k) =  dtauy(i,k) + a(i,k)
      b(i,k) =  dtaux(i,k) + b(i,k)
      dusfc(i)=dusfc(i)+dtaux(i,k)*del(i,k)
      dvsfc(i)=dvsfc(i)+dtauy(i,k)*del(i,k)
 2000 continue
      do 2002 i=1,idimt
      dusfc(i)=-1.e3/g*rcs*pstar(i)*dusfc(i)
      dvsfc(i)=-1.e3/g*rcs*pstar(i)*dvsfc(i)
 2002 continue
c
c-----diagnostic flag .ne.0 on, otherwise off, num(715)=output unit (=6)
c
      return
      end
