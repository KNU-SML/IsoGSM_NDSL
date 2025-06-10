#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine swflux99(ipts,tau,ssc,g0,ff,csm,zth,alb,ald,
     1                  daytm,upflux,dwflux,dwsfcb,dwsfcd)
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
c    ff:  the effective forward scattering factor
c    csm: secant of the zenith angle
c    zth: cosine of the zenith angle
c    alb: surface albedo for direct radiation
c    ald: surface albedo for diffused radiation
c    daytm: daytime flag
c
c  outputs:
c    upflux: upward fluxes
c    dwflux: downward fluxes
c    dwsfcb: downward surface flux direct component
c    dwsfcd: downward surface flux diffused component
c********************************************************************
c
#include <paramodel.h>
#include <rdparm99.h>
c --- input
                       d i m e n s i o n
     1  tau(imbx,l), ssc(imbx,l), g0(imbx,l), ff(imbx,l)
     2, csm(imax),   zth(imax),   alb(imax),  ald(imax)
      logical daytm(imax)
c --- output
                       d i m e n s i o n
     1  upflux(imbx,lp1),dwflux(imbx,lp1),dwsfcb(imax),dwsfcd(imax)
c --- temporary
                       d i m e n s i o n
     1  ttb(imbx,lp1),tdn(imbx,lp1),rup(imbx,lp1), tt (imbx,lp1,2)
     2, rfu(imbx,lp1),rfd(imbx,lp1),tb (imbx,lp1), rr (imbx,lp1,2)
c
c===> ... diffuse incident radiation is approximated by beam radiation
c         with an incident angle of 53 degrees. cos(53) = 0.602
      zthd = 0.602
      csmd = 1.0 / zthd
      epsln = 1.0e-30
c
c===> ... delta-eddington scaling of single scattering albedo,
c         optical thickness, and asymmetry factor, k & h eqs(27-29)
c
      do k=1,l
      do i=1,ipts
        if (daytm(i)) then
c===> ... delta-eddington scaling of single scattering albedo,
c         optical thickness, and asymmetry factor, k & h eqs(27-29)
          aa   = 1.0e0 - ff(i,k)*ssc(i,k)
          taup = tau(i,k) * aa
          sscp = ssc(i,k) * (1.0e0 - ff(i,k)) / aa
          gp   = (g0(i,k) - ff(i,k)) / (1.0e0 - ff(i,k))
c
          oms1 = 1.0e0 - sscp
          ogs1 = 1.0e0 - sscp*gp
          tlam = 3.0e0 * oms1*ogs1
          slam = sqrt(tlam)

          u1   = 1.50e0 * ogs1 / slam
          u1p1 = u1 + 1.0e0
          u1m1 = u1 - 1.0e0
          e1   = exp(max(-taup*slam, -30.0))
          u1e  = u1 * e1
          u1epe= u1e + e1
          u1eme= u1e - e1
          den  = 1.0e0 / ((u1p1 + u1eme)*(u1p1 - u1eme))
          rf1  = (u1p1 + u1epe) * (u1m1 - u1eme) * den
          tf1  = 4.0e0 * u1e * den

c===> ... compute layer transmissions and reflections
c         (i,k,j) j=1,2 for layer k illuminated by diffuse and
c                       direct incoming radiation
c         rr   :  layer reflection
c         tt   :  layer total transmission
c         tb   :  layer direc transmission
c
c       diffuse radiation
c       -----------------
          zzth = zthd
          zz   = zzth * zzth
          den1 = 1.0 - tlam*zz
c===> ... safety check
          if (abs(den1) .lt. 1.0e-8) then
            zzth = zzth + 0.001
            zz   = zzth * zzth
            den1 = 1.0 - tlam*zz
          end if
          den1 = sscp / den1
c
          gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den1
          alfa = 0.75e0 * zthd * (gp + ogs1) * den1
          amg  = alfa - gama
          apg  = alfa + gama
c
          tb(i,k)  = exp( -min(30.0, taup*csmd) )
          za       = amg * tb(i,k)
          rr(i,k,1)= za*tf1 + apg*rf1 - amg
          tt(i,k,1)= za*rf1 + apg*tf1 + (1.0e0-apg)*tb(i,k)
c
c       direct radiation
c       -----------------
          zzth = zth(i)
          zz   = zzth * zzth
          den1 = 1.0 - tlam*zz
c===> ... safety check
          if (abs(den1) .lt. 1.0e-8) then
            zzth = zzth + 0.001
            zz   = zzth * zzth
            den1 = 1.0 - tlam*zz
          end if
          den1   = sscp / den1
c
          gama = 0.50e0 * (1.0e0 + 3.0e0*gp*oms1*zz) * den1
          alfa = 0.75e0 * zth(i) * (gp + ogs1) * den1
          amg  = alfa - gama
          apg  = alfa + gama
c
          tb(i,k)  = exp( -min(30.0, taup*csm(i)) )
          za       = amg * tb(i,k)
          rr(i,k,2)= za*tf1 + apg*rf1 - amg
          tt(i,k,2)= za*rf1 + apg*tf1 + (1.0e0-apg)*tb(i,k)
c
          tb(i,k)   = max(0.0e0, tb(i,k))
          rr(i,k,2) = max(0.0e0, rr(i,k,2))
          tt(i,k,2) = max(0.0e0, tt(i,k,2))
          rr(i,k,1) = max(0.0e0, rr(i,k,1))
          tt(i,k,1) = max(0.0e0, tt(i,k,1))
c --- night time condition
        else
          tb(i,k)   = 0.0e0
          rr(i,k,1) = 0.0e0
          rr(i,k,2) = 0.0e0
          tt(i,k,1) = 0.0e0
          tt(i,k,2) = 0.0e0
        end if
      end do
      end do
c
c --- at the surface
c
      do i=1,ipts
        tb(i,lp1)   = 0.0e0
        rr(i,lp1,2) = alb(i)
        tt(i,lp1,2) = 0.0e0
        rr(i,lp1,1) = ald(i)
        tt(i,lp1,1) = 0.0e0
      end do
c
      do i=1,ipts
        ttb(i,1) = tb(i,1)
        tdn(i,1) = tt(i,1,2)
        rfd(i,1) = rr(i,1,1)
      end do
c
c===> ... layers added downward starting from top
c
      do k=2,lp1
      do i=1,ipts
        if (daytm(i)) then
          den = tt(i,k,1) / (1.0e0 - rfd(i,k-1) * rr(i,k,1))
          ttb(i,k) = ttb(i,k-1) * tb(i,k)
          if (ttb(i,k) .lt. epsln) ttb(i,k) = 0.0
          tdn(i,k) = ttb(i,k-1)*tt(i,k,2)+(tdn(i,k-1)-ttb(i,k-1)
     1             + ttb(i,k-1)*rr(i,k,2)*rfd(i,k-1)) * den
          rfd(i,k) = rr(i,k,1) + tt(i,k,1)*rfd(i,k-1) * den
        end if
      end do
      end do
c
c===> ... layers added upward starting from surface
c
      do i=1,ipts
        rfu(i,lp1) = rr(i,lp1,1)
        rup(i,lp1) = rr(i,lp1,2)
      end do
      do k=l,1,-1
        kp1 = k + 1
        do i=1,ipts
          if (daytm(i)) then
            den = tt(i,k,1) / (1.0e0 - rfu(i,kp1) * rr(i,k,1))
            rup(i,k) = rr(i,k,2) + ((tt(i,k,2)-tb(i,k))*rfu(i,kp1)
     1               + tb(i,k)*rup(i,kp1)) * den
            rfu(i,k) = rr(i,k,1) + tt(i,k,1)*rfu(i,kp1) * den
          end if
        end do
      end do
c
c===> ... find upward and downward fluxes
c
      do i=1,ipts
        if (daytm(i)) then
          upflux(i,1) = rup(i,1)
          dwflux(i,1) = 1.0e0
        else
          upflux(i,1) = 0.0e0
          dwflux(i,1) = 0.0e0
        end if
      end do
      do k=2,lp1
        km1 = k - 1
        do i=1,ipts
          if (daytm(i)) then
            den = 1.0e0 / (1.0e0 - rfd(i,km1)*rfu(i,k))
            aa  = ttb(i,km1) * rup(i,k)
            bb  = tdn(i,km1) - ttb(i,km1)
            upflux(i,k) = (aa + bb*rfu(i,k)) * den
            dwflux(i,k) = ttb(i,km1) + (aa*rfd(i,km1) + bb) * den
          else
            upflux(i,k) = 0.0e0
            dwflux(i,k) = 0.0e0
          end if
        end do
      end do
c
c===> ... surface downward fluxes
c
      do i=1,ipts
        dwsfcb(i) = ttb(i,l)
        dwsfcd(i) = dwflux(i,lp1)-dwsfcb(i)
      end do
c
      return
      end
