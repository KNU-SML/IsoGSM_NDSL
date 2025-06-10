#define PI pi
#define RERTH rerth
      subroutine ll2xyr(io,cproj,corient,ctruth,ccotru,clat,clon,cx,cy)
c
#include <paramodel.h>
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.,rerth=rerth_)
c
      real clat(io),clon(io),cx(io),cy(io)
c
c input are all degree and output to global x y, not domain x y
c
c if proj= 0  do mercater projection
c if proj= 1  do north polar projection
c if proj=-1  do south polar projection
c if proj= 2  do north lambert projection
c if proj=-2  do south lambert projection
c if proj= 5  do lambert azimuthal equal area projection
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
          blat = clat(i) * rad
          blon = clon(i) * rad
          rsoa2 = tan( (hfpi-blat*nproj)*0.5 )
          x2py2 = ( rsoa2 * a2 ) ** 2.0
          blon = mod(blon,twopi)
          if(blon.lt.0. e 0) blon = twopi + blon
          rlon = nproj * (blon - dlamda0)
          rlon = amod(rlon,twopi)
          if( rlon.lt.0. ) rlon=twopi+rlon
          yox = tan(rlon)
          x = sqrt( x2py2/(1.+yox*yox) )
          y = sqrt( x2py2 - x*x )
          if( rlon.gt.hfpi .and. rlon.lt. pi+hfpi ) x = -x
          if( rlon.gt.pi .and. rlon.lt. twopi ) y = -y
          cx(i) = x
          cy(i) = y
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
          blat = clat(i) * rad
          blon = clon(i) * rad

c regional mercator projection that crosses 0E.
c hyodae/kana 6/3/2008

       if(io.gt.1) then
       	  if(clon(1).lt.clon(io)) then
            blon = mod(blon,twopi)
          else
c if clon(1) is greater than clon(end)
c and if clon(1)>0
c then the domain crosses the dateline (eg. 270E-20E)
c if so, add twopi to clon(i) EAST of 0E (e.g. 270E-380E)
         if(clon(i).ge.0.and.clon(i).le.clon(io)) then
          blon=blon+twopi
         endif
        endif
      else
          blon = mod(blon,twopi)
      endif

c          blon = mod(blon,twopi)

          if(blon.lt.0. e 0) blon = twopi + blon
          x=a2*(blon-cenlon)
          y=a2*log(tan(blat/2.0+qtpi))
          cx(i) = x
          cy(i) = y
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
          blat = clat(i) * rad
          blon = clon(i) * rad
          blon = mod(blon,twopi)
          if(blon.lt.0. e 0) blon = twopi + blon
          r=r00*(tan(qtpi-is*blat/2))**cone
          x=    r*sin(cone*(blon-cenlon))
          y=-is*r*cos(cone*(blon-cenlon))
          cx(i) = x
          cy(i) = y
        enddo
c
c  lambert szinuthal equal area
c
      else if( nproj.eq.5) then
        alamd0  = corient * rad
        phi0    = ctruth * rad
        cosphi0 = cos(phi0)
        sinphi0 = sin(phi0)
        do i=1,io
          blat = clat(i) * rad
          blon = (clon(i) - corient) * rad
          akp=1.+sinphi0*sin(blat)+cosphi0*cos(blat)*cos(blon)
          akp=sqrt(2./akp)
          cx(i)=akp*cos(blat)*sin(blon)*rerth_
          cy(i)=akp*(cosphi0*sin(blat)-sinphi0*cos(blat)*cos(blon))*
     2          rerth_
        enddo
      endif
c
      return
      end
