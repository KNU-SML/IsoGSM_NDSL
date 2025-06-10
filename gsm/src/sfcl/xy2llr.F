#define PI pi
#define RERTH rerth
      subroutine xy2llr(io,cproj,corient,ctruth,ccotru,cx,cy,clat,clon)
c
#include <paramodel.h>
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.,rerth=rerth_)
c
      real cx(io),cy(io)
      real clat(io),clon(io)
c
c input are all degree and output to global x y, not regional x y
c
c if proj=0  do mercator projection
c if proj=1  do north polar projection
c if proj=-1 do south polar projection
c if proj=2  do north lambert projection
c if proj=-2 do south lambert projection
c if proj=5  do lambert azimuthal equal area projection
c
c polar projection
c
      nproj = cproj
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
        truth  = ctruth * rad
        truth  = nproj * truth
        orient  = corient * rad
        cenlon = mod(orient,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        dlamda0 = cenlon + hfpi
        a2 =  rerth * ( 1.0 + sin(truth) )
        radlat = 90. * rad
        radlon = 0.0 * rad - dlamda0
        radlat = nproj * radlat
        radlon = nproj * radlon
        do i=1,io
          x = cx(i)
          y = cy(i)
          if( x.gt.0.0 ) then
            blon = atan(y/x)
          else if( x.lt.0.0 ) then
            blon = pi + atan(y/x)
          else
            blon = hfpi
            if( y.lt.0.0 ) blon = blon * 3.0
          endif
ckeic          blon = blon + dlamda0
          blon = blon + nproj*dlamda0
          blon = mod(blon,twopi)
          blon = nproj * blon
          rsoa2 = sqrt( x*x + y*y )/a2
ckeic          blat = hfpi - 2. * atan(rsoa2)
          blat = nproj * ( hfpi - 2. * atan(rsoa2) )
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
c
c  mercator
c
      else if ( nproj.eq.0 ) then
        truth  = ctruth * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        a2 =  rerth * cos( truth )
        dlamda0 = 0.0
        do i=1,io
          x = cx(i)
          y = cy(i)
          blon = x / a2 + cenlon
          blon = mod(blon,twopi)
          if(blon.lt.0. e 0) blon = twopi + blon
          blat = 2.*(atan(exp(y/a2))-qtpi)
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
c
c  lambert
c
      else if( nproj.eq.2 .or. nproj.eq.-2 ) then
        is=1
        if( nproj.lt.0 ) is=-1
        truth  = ctruth * rad
        cotru  = ccotru * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        if( ctruth.eq.ccotru ) then
          cone= cos (hfpi-is*truth)
        else
          cone=(log(cos(truth))-log(cos(cotru)))/
     1         (log(tan(qtpi-is*truth/2))-log(tan(qtpi-is*cotru/2)))
        endif
        dlamda0 = 0.0
        r00=rerth/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
        do i=1,io
          r = sqrt(cx(i)*cx(i) + cy(i)*cy(i) )
          blon = cenlon + asin(cx(i)/r) / cone
          blon = mod(blon+twopi,twopi)
          blat = hfpi - 2 * is * atan ( (r/r00)**(1./cone) )
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
c
c  lambert azimuthal eual area
c    corient is the central longitude (e.g. -100. over us)
c    ctruth  is the latitude of origin or standard parallel (e.g. 45. over us)
c
      else if(nproj.eq.5 ) then
        alamd0  = corient * rad
        phi0    = ctruth * rad
        cosphi0 = cos(phi0)
        sinphi0 = sin(phi0)
        do i=1,io
          r=sqrt(cx(i)*cx(i) + cy(i)*cy(i))
          if(0.5*r/rerth_.gt.1) then
             print *,'error in Lambert Azimuthal'
             print *,'0.5*r/rerth_=',0.5*r/rerth_
             call abort
          endif
          c=2.*asin(0.5*r/rerth_)
          sinc=sin(c)
          cosc=cos(c)
          clat(i)=asin(cosc*sinphi0+cy(i)*sinc*cosphi0/r)
          clon(i)=alamd0+atan(cx(i)*sinc/
     1              (r*cosphi0*cosc-cy(i)*sinphi0*sinc))
          clat(i) = clat(i) / rad
          clon(i) = clon(i) / rad
        enddo
      endif
c
      return
      end
