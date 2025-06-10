#define IGRD igrd
#define IGRD1 igrd1
#define IGRD12 igrd12
#define JGRD jgrd
#define JGRD1 jgrd1
#define JGRD12 jgrd12
#define LNGRD lngrd
#define PI pi
#define RERTH rerth
#define NPES npes
#define G2G
      subroutine setgrd(flat,flon,delx,dely,dlamda0)
#include <paramodel.h>
#ifndef G2G
      dimension flat(lngrd_),flon(lngrd_)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  setgrd
c   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
c   prgmmr:  yifeng cui, mp print       org: sdsc       date: 02-09-11
c
c abstract:  setup grid latitude, longitude, and map factor etc for
c            regional grid.
c
c program history log:
c
c usage:    call setgrd(flat,flon,delx,dely,dlamda0)
c   output argument list:
c     flat	- latitude
c     flon	- logitude
c     delx	- grid spacing in x direction
c     dely	- grid spacing in y direction
c     dlamda0	- angle between 0 longitude and x positive axies
c
c   input argument list:
c
c   input files: none
c
c   output files: none
c
c   subprograms called: none
c
c   remark: none
c
c attributes:
c   language: fortran 77.
c   machine:  cray ymp.
c
c$$$

#include <rscomloc.h>
c#ifdef MP
c#include <comrmpi.h>
c#endif
c
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.)
c
#ifdef DBG
c#ifdef MP
c      if (mype.eq.master) then
c#endif
      print 1234
 1234 format(' ==== in routine setgrd === ')
c#ifdef MP
c      endif
c#endif
#endif
c
      delx = rdelx
      dely = rdely
c
c --------- setup regional lat/lon and map factor -----
c
c if proj=0  do mercater projection
c if proj=1  do north polar projection
c if proj=-1 do south polar projection
c
      nproj = rproj
c
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
c ++++++++++++++++++++++++++++++++++++++
c polar projection
c ++++++++++++++++++++++++++++++++++++++
      truth  = rtruth * rad
      truth  = nproj * truth
      orient  = rorient * rad
      dlamda0 = orient + hfpi
      a2 =  rerth_ * ( 1.0 + sin(truth) )
      radlat = rcenlat * rad
      radlon = rcenlon * rad - dlamda0
      radlat = nproj * radlat
      radlon = nproj * radlon
      yyy = a2 * cos( radlat )/(1. + sin( radlat ) )
      cenlat = rcenlat
      if( abs(cenlat) .eq. 90. ) then yyy = 0.0
      y00 = yyy * sin( radlon ) - ( rbtmgrd -1.) * dely
      x00 = yyy * cos( radlon ) - ( rlftgrd -1.) * delx
#ifdef DBG
c#ifdef MP
c      if (mype.eq.master) then
c#endif
      print *,' delx x00 y00 ',delx,x00,y00
c#ifdef MP
c      endif
c#endif
#endif
c
c =========
c           lat loop
      do 100 j = 1,jgrd12_
      lats = j
      ijlats = (lats-1)*igrd12_
      ys = y00 + (lats-1)*dely
      latn = jgrd1_+1-j
      ijlatn = ijlats+igrd1_
      yn = y00 + (latn-1)*dely
c
      do 100 i=1,igrd1_
      x = x00 + (i-1)*delx
      if( x .gt. 0. e 0 ) then
         flons = atan(ys/x)
         flonn = atan(yn/x)
      else if ( x .lt. 0. e 0 ) then
         flons = pi + atan(ys/x)
         flonn = pi + atan(yn/x)
      else
         flons = hfpi
         if( ys .lt. 0. e 0 ) flons = flons * 3.0
         flonn = hfpi
         if( yn .lt. 0. e 0 ) flonn = flonn * 3.0
      endif
      flons = nproj * flons + dlamda0
      flons = mod(flons,twopi)
      if(flons.lt.0. e 0) flons = twopi + flons
      flonn = nproj * flonn + dlamda0
      flonn = mod(flonn,twopi)
      if(flonn.lt.0. e 0) flonn = twopi + flonn
c
      rsoa2 = sqrt( x*x + ys*ys )/a2
      flats = hfpi - 2.0 e 0 * atan(rsoa2)
      flat(ijlats+i) = nproj * flats
      flon(ijlats+i) = flons
c
      rsoa2 = sqrt( x*x + yn*yn )/a2
      flatn = hfpi - 2.0 e 0 * atan(rsoa2)
      flat(ijlatn+i) = nproj * flatn
      flon(ijlatn+i) = flonn
c
 100  continue
c
      else if ( nproj.eq.0 ) then
c
c ++++++++++++++++++++++++++++
c do mercater
c ++++++++++++++++++++++++++++
      truth  = rtruth * rad
      cenlat = rcenlat * rad
      cenlon = rcenlon * rad
      cenlon = mod(cenlon,twopi)
      if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
      a2 =  rerth_ * cos( truth )
      x0 = 0.0
      y0 = a2 * log( abs( tan( qtpi + 0.5 * cenlat ) ) )
      x00 = - ( rlftgrd - 1.0 ) * delx                           
      y00 = - ( rbtmgrd - 1.0 ) * dely
      dlamda0 = 0.0
c
      do 200 j = 1,jgrd12_
      lats = j
      ijlats = (lats-1)*igrd12_
      ys = y00 + (lats-1)*dely + y0
      latn = jgrd1_+1-j
      ijlatn = ijlats+igrd1_
      yn = y00 + (latn-1)*dely + y0
c
       do 200 i=1,igrd1_
         x = x00 + (i-1)*delx + x0
         flons = x / a2 + cenlon
         flons = mod(flons,twopi)
         if(flons.lt.0. e 0) flons = twopi + flons
         flonn = flons
c
         flats = 2.0 *( atan( exp( ys/a2 ) ) - qtpi )
         flat(ijlats+i) = flats
         flon(ijlats+i) = flons
c
         flatn = 2.0 *( atan( exp( yn/a2 ) ) - qtpi )
         flat(ijlatn+i) = flatn
         flon(ijlatn+i) = flonn
c
 200  continue
c
      endif
c
c
#endif
      return
      end
