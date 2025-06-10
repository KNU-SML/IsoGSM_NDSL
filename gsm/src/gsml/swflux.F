#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine swflux(ipts,tau,ssc,g0,csm,alb,ald,daytm,
     1                  isbd,upflux,dwflux,dwsfcb,dwsfcd)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c********************************************************************
c  uses the delta-eddington approximation to compute the bulk
c  scattering properties of a single layer coded following
c  coakley et al.  (jas, 1982)
c
c  inputs:
c    tau: the effective optical thickness
c    ssc: the effective single scattering albedo
c    g0:  the effective asymmetry factor
c    csm: the effective secant of the zenith angle
c    alb: surface albedo for direct radiation
c    ald: surface albedo for diffused radiation
c    daytm: daytime flag
c    isbd: =1 for uv+vis spectral bands
c          =2 for nir spectral bands
c
c  outputs:
c    upflux: upward fluxes
c    dwflux: downward fluxes
c    dwsfcb: downward surface flux direct component
c    dwsfcd: downward surface flux diffused component
c********************************************************************
c
#include <paramodel.h>
#include <rdparm.h>
c --- input
                       d i m e n s i o n
     1  tau(imbx,l), ssc(imbx,l), g0(imbx,l)
     2, csm(imax),   alb(imax),   ald(imax)
      logical daytm(imax)
c --- output
                       d i m e n s i o n
     1  upflux(imbx,lp1),dwflux(imbx,lp1),dwsfcb(imax),dwsfcd(imax)
c --- temporary
                       d i m e n s i o n
     1  tb (imbx,lp1),r0 (imbx,lp1),t0 (imbx,lp1), rf (imbx,lp1)
     2, ttb(imbx,lp1),tdn(imbx,lp1),rup(imbx,lp1), tf (imbx,lp1)
     3, tfd(imbx,lp1),rfu(imbx,lp1),rfd(imbx,lp1), zth(imax)
c
      do 10 i=1,ipts
        tb(i,lp1) = 0.0e0
        r0(i,lp1) = alb(i)
        t0(i,lp1) = 0.0e0
        rf(i,lp1) = ald(i)
        tf(i,lp1) = 0.0e0
        zth(i)    = 1.0e0 / csm(i)
  10  continue
c
      do 20 k=1,l
      do 20 i=1,ipts
       if (daytm(i)) then
c===> ... delta-eddington scaling of single scattering albedo,
c         optical thickness, and asymmetry factor, k & h eqs(27-29)
        ff   = g0(i,k) * g0(i,k)
        aa   = 1.0e0 - ff*ssc(i,k)
        taup = tau(i,k) * aa
        sscp = ssc(i,k) * (1.0e0 - ff) / aa
        gp   = g0(i,k) / (1.0e0 + g0(i,k))
c
        oms1 = 1.0e0 - sscp
        ogs1 = 1.0e0 - sscp*gp
        tlam = 3.0e0 * oms1*ogs1
        slam = sqrt(tlam)
        zz   = zth(i) * zth(i)
        den1 = 1.0e0 - tlam*zz
c===> ... safety check
        den  = sscp / sign(amax1(1.0e-20, abs(den1)), den1)
c
        gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den
        alfa = 0.75e0 * zth(i) * (gp + ogs1) * den
        u1   = 1.50e0 * ogs1 / slam
        up1  = u1 + 1.0e0
        um1  = u1 - 1.0e0
        amg  = alfa - gama
        apg  = alfa + gama
c
c===> ... compute layer transmissions and reflections
c         r0   :  layer reflection of the direct beam
c         t0   :  layer diffuse+direct transmission of direct beam
c         rf   :  layer reflection of the diffused radiation
c         tf   :  layer transmission of the diffused radiation
c         tb   :  layer direc transmission of the direct beam
c
        e1  = exp( -taup*slam )
        ue  = u1 * e1
        uepe= ue + e1
        ueme= ue - e1
        den = 1.0e0 / ((up1 + ueme)*(up1 - ueme))
c
        arg = amin1(30.0e0, taup*csm(i))
        tb1  = exp(-arg)
        rf1  = (up1 + uepe) * (um1 - ueme) * den
        tf1  = 4.0e0 * ue * den
        za   = amg * tb1
        r01  = za * tf1 + apg * rf1 - amg
        t01  = za * rf1 + apg * tf1 - (apg - 1.0e0)*tb1
c
        tb(i,k) = amax1(0.0e0, tb1)
        r0(i,k) = amax1(0.0e0, r01)
        t0(i,k) = amax1(0.0e0, t01)
        rf(i,k) = amax1(0.0e0, rf1)
        tf(i,k) = amax1(0.0e0, tf1)
       else
        tb(i,k) = 0.0e0
        r0(i,k) = 0.0e0
        t0(i,k) = 0.0e0
        rf(i,k) = 0.0e0
        tf(i,k) = 0.0e0
       end if
  20  continue
c
      if (isbd .eq. 2) then
        do 40 k=1,l
        do 40 i=1,ipts
          if (daytm(i).and.ssc(i,k).le.0.0001e0) then
            tb(i,k) = exp(-tau(i,k)*csm(i))
            t0(i,k) = tb(i,k)
            r0(i,k) = 0.0e0
            tf(i,k) = exp(-1.66e0*tau(i,k))
            rf(i,k) = 0.0e0
          end if
  40    continue
      end if
c
      do 60 i=1,ipts
        tdn(i,1) = t0(i,1)
        rfd(i,1) = rf(i,1)
        tfd(i,1) = tf(i,1)
        ttb(i,1) = tb(i,1)
        ttb(i,l) = 0.0e0
  60  continue
c
c===> ... layers added downward starting from top
c
      do 80 k=2,lp1
      do 80 i=1,ipts
       if (daytm(i)) then
        den = tf(i,k) / (1.0e0 - rfd(i,k-1) * rf(i,k))
        tdn(i,k) = ttb(i,k-1)*t0(i,k) + (tdn(i,k-1)-ttb(i,k-1)
     1           + ttb(i,k-1)*r0(i,k)*rfd(i,k-1)) * den
        rfd(i,k) = rf(i,k) + tf(i,k)*rfd(i,k-1) * den
        tfd(i,k) = tfd(i,k-1) * den
        if(abs(tfd(i,k)).lt.1.e-20) tfd(i,k)=0.
        ttb(i,k) = ttb(i,k-1) * tb(i,k)
        if(abs(ttb(i,k)).lt.1.e-20) ttb(i,k)=0.
       end if
  80  continue
c
c===> ... layers added upward starting from surface
c
      do 100 i=1,ipts
        rfu(i,lp1) = rf(i,lp1)
        rup(i,lp1) = r0(i,lp1)
 100  continue
      do 120 k=l,1,-1
      do 120 i=1,ipts
       if (daytm(i)) then
        den = tf(i,k) / (1.0e0 - rfu(i,k+1) * rf(i,k))
        rup(i,k) = r0(i,k) + ((t0(i,k)-tb(i,k))*rfu(i,k+1)
     1           + tb(i,k)*rup(i,k+1)) * den
        rfu(i,k) = rf(i,k) + tf(i,k)*rfu(i,k+1) * den
       end if
 120  continue
c
c===> ... find upward and downward fluxes
c
      do 160 i=1,ipts
       if (daytm(i)) then
        upflux(i,1) = rup(i,1)
        dwflux(i,1) = 1.0e0
       else
        upflux(i,1) = 0.0e0
        dwflux(i,1) = 0.0e0
       end if
 160  continue
      do 180 k=2,lp1
      do 180 i=1,ipts
       if (daytm(i)) then
        den = 1.0e0 / (1.0e0 - rfd(i,k-1)*rfu(i,k))
        upflux(i,k) = (ttb(i,k-1)*rup(i,k) +
     1                (tdn(i,k-1)-ttb(i,k-1))*rfu(i,k)) * den
        dwflux(i,k) = ttb(i,k-1) + ((tdn(i,k-1)-ttb(i,k-1))
     1              + ttb(i,k-1)*rup(i,k)*rfd(i,k-1)) * den
       else
        upflux(i,k) = 0.0e0
        dwflux(i,k) = 0.0e0
       end if
 180  continue
c
c===> ... surface downward fluxes
c
      do 200 i=1,ipts
        dwsfcb(i) = ttb(i,l)
        dwsfcd(i) = dwflux(i,lp1)-dwsfcb(i)
 200  continue
c
      return
      end
