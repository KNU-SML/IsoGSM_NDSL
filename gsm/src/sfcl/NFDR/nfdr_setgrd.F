#define PI pi
#define RERTH rerth
#define NPES npes
      subroutine nfdr_setgrd(rlat,rlon,delx,dely)
c 
      implicit none
c
c
      integer idim,jdim
      parameter(idim=184,jdim=116)
      real rproj,rtruth,rorient,rdelx,rdely
      parameter(rproj=5.,rtruth=45.,rorient=-100.,
     1          rdelx=25000.,rdely=25000.)
      real rlftgrd,rbtmgrd
      parameter(rlftgrd=-82.,rbtmgrd=-85.)
c
      real rlat(idim,jdim),rlon(idim,jdim),delx,dely
c
#include <paramodel.h>
c
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
c
      real pi,twopi,hfpi,qtpi,rad
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.)
c
      real alamd0,phi0,cosphi0,sinphi0,x,y,r,c,sinc,cosc
      integer i,j
      real x00,y00
c
c      delx = rdelx
c      dely = rdely
c
c --------- setup regional lat/lon and map factor -----
c
c this is for Lambert azimuthal equal area projection only
c
      x00 = rlftgrd * delx 
      y00 = rbtmgrd * dely
      alamd0  = rorient * rad
      phi0    = rtruth * rad
      cosphi0 = cos(phi0)
      sinphi0 = sin(phi0)
      do j = 1,jdim
        y = y00 + (j-1)*dely
        do i=1,idim
          x = x00 + (i-1)*delx
          r=sqrt(x*x + y*y)
          if(0.5*r/rerth_.gt.1.) then
            print *,'error in nfdr_setgrd'
            print *,'0.5*r/rerth_=',0.5*r/rerth_
            call abort
          endif
          c=2.*asin(0.5*r/rerth_)
          sinc=sin(c)
          cosc=cos(c)
          rlat(i,j)=asin(cosc*sinphi0+y*sinc*cosphi0/r)
          rlon(i,j)=alamd0+atan(x*sinc/
     1              (r*cosphi0*cosc-y*sinphi0*sinc))
          rlat(i,j) = rlat(i,j) / rad
          rlon(i,j) = rlon(i,j) / rad
        enddo
      enddo
c
c     print *,'rlat=',rlat
c     print *,'rlon=',rlon
c
      return
      end
