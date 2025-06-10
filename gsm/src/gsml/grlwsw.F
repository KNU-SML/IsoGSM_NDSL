#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
#define PI pi
      subroutine grlwsw(sigl,nfile)
#include <paramodel.h>
#include <rdparm.h>
      parameter (lngth=37*l)
      parameter (nl=81,nlp1=nl+1,nlgth=37*nl)
      common /rdfsav/ degrad,hsigma,daysec,rco2
c   **************************************************************
c   seasonal climatologies of o3 (obtained on user vertical coord)
c   defined as 5 deg lat means n.p.->s.p.
c   **************************************************************
      common /savmem/
c ---   ...winter....  ...spring....  ...summer....  ....fall.....
     1   dduo3n(37,l), ddo3n2(37,l), ddo3n3(37,l), ddo3n4(37,l)
      dimension rad1(lngth), rad2(lngth), rad3(lngth), rad4(lngth)
      equivalence (rad1(1),dduo3n(1,1)),(rad2(1),ddo3n2(1,1))
      equivalence (rad3(1),ddo3n3(1,1)),(rad4(1),ddo3n4(1,1))
c ---  seasonal climatologies of o3 on the detailed gfdl coordinate...
c      interpolation to each point pressure profile done in ozon2d
c      (see gloopr)
      common /seaso3/
c ---   ...winter....  ...spring....  ...summer....  ....fall.....
     1   xduo3n(37,nl), xdo3n2(37,nl), xdo3n3(37,nl), xdo3n4(37,nl)
     2,  prgfdl(nl)
      dimension xrad1(nlgth),xrad2(nlgth),xrad3(nlgth),xrad4(nlgth)
      equivalence (xrad1(1),xduo3n(1,1)),(xrad2(1),xdo3n2(1,1))
      equivalence (xrad3(1),xdo3n3(1,1)),(xrad4(1),xdo3n4(1,1))
c
      dimension sigl(l),pstd(nlp1)
c   ********************************************************
c   *      one time computation of necessary quantities    *
c   ********************************************************
c===> ... initialize arrays,get constants,etc...
      degrad=180.0e0/pi_
      hsigma=5.673e-5
      daysec=1.1574e-5
c===> ... atmosperic carbon dioxide concentration is now read by conrad
c         but it defaults to 348 ppm for backward compatibility.
#ifndef CO2PPM
      rco2=3.48e-4
#else
      rco2=CO2PPMe-6
#endif
      call hconst
c===> ... interpolate climo o3 to the current vertical coordinate
c         need layer sigma, get from psfc and layer p for i=1
      call o3int(dduo3n,sigl)
c===> ... compute detailed o3 profile from the original gfdl pressures
c         where output from o3int (pstd) is top down in mb*1.e3
c         and psfc=1013.25 mb    ......k.a.c. dec94
      call o3intn(xduo3n,pstd)
#ifndef CO2PPM
      call conrad(nfile,rco2)
#endif
c===> ... average climatological valus of o3 from 5 deg lat means,
c         so that time and space interpolation will work
      do 5 i=1,lngth
        avg=.25e0*(rad1(i)+rad2(i)+rad3(i)+rad4(i))
        a1=.5e0*(rad2(i)-rad4(i))
        b1=.5e0*(rad1(i)-rad3(i))
        b2=.25e0*((rad1(i)+rad3(i))-(rad2(i)+rad4(i)))
        rad1(i)=avg
        rad2(i)=a1
        rad3(i)=b1
        rad4(i)=b2
    5 continue
c===> ... average climatological valus of o3 from 5 deg lat means,
c         so that time and space interpolation will work
c        (see subprogram ozon2d)
      do 10 i=1,nlgth
        avg=.25e0*(xrad1(i)+xrad2(i)+xrad3(i)+xrad4(i))
        a1=.5e0*(xrad2(i)-xrad4(i))
        b1=.5e0*(xrad1(i)-xrad3(i))
        b2=.25e0*((xrad1(i)+xrad3(i))-(xrad2(i)+xrad4(i)))
        xrad1(i)=avg
        xrad2(i)=a1
        xrad3(i)=b1
        xrad4(i)=b2
   10 continue
c===> ... get gfdl pressure in cb (flip vertical coordinate)
      do 20 n=1,nl
       prgfdl(n) = pstd(nl+1-n)*1.e-4
   20 continue
c
      return
      end
