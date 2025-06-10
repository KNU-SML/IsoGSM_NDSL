#define JCAP jcap
#define LATG latg
#define NPES npes
      subroutine o3sbuv(fhour,idate,ko3,o3out,pstr,jerr)
c     ********************************************************
c     *  computes o3 climo from 12 month dataset, linearly   *
c     *   interpolated to day,mon of the fcst.  then create  *
c     *   a 5 deg array from the 10 deg climatology...for    *
c     *   ease when doing a latitudinal interpolation        *
c     *  thanks to s moorthi for new o3 climo...kac  dec 1996*
c     * input:                                               *
c     *   idate=nmc date-time                                *
c     *   fhour=forecast hour                                *
c     *   ko3=unit number of o3 climatology                  *
c     * output :                                             *
c     *   o3out=5-deg o3 climo for forecast date(np->spole)  *
c     *   pstr=pressure (mb) for the climo lyrs (k=1 is top) *
c     *   jerr=0 if o3 file exists, =1 if not (gfdl=default) *
c     ********************************************************
c     * yifeng cui (sdsc) added rsm mp    date:02-09-05      *
c     ********************************************************
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#else
#include <commpi.h>
#endif
#endif
c geos ozone data
      parameter (jmr=18, blte=-85.0, dlte=10.0)
      parameter (loz=17,jmout=37)
c
      integer   imon(12),days(12),idate(4)
      dimension  ilat(jmr,12)
      dimension o3r(jmr,loz,12),o3tmp(jmr,loz),o3tmp2(jmr,loz)
      dimension o3out(jmout,loz),pstr(loz)
      data  days/31,28,31,30,31,30,31,31,30,31,30,31/
c
      ida=idate(3)
      imo=idate(2)
c
c   find current day and month, initial values in ida,imo!
c       will not worry about leap year, since it will take a
c       120-year (what?) forecast to be off by 1 month.  if this
c       is deemed a problem, need to redo this calculation.
c
      if (fhour.ge.24.) then
c
c  number of days into the forecast
c
       numdyz=fhour/24.0 + 0.01
c
c  get day-of-year, remember climate runs are for years
c
       imo1=imo-1
       jday = ida
       if (imo1.gt.0) then
        jday=0
        do 7 ken=1,imo1
         jday=jday+days(ken)
    7   continue
        jday=jday+ida
       end if
       nmdtot = jday+numdyz
       ndayr = mod(nmdtot,365)
       if (ndayr.eq.0) ndayr=365
c
c   now get month from day-of-year
c
       mday=0
       do 8 ken=1,11
        mday=mday+days(ken)
        imo=ken
        if (ndayr.le.mday) then
         ida=ndayr-(mday-days(imo))
         go to 9
        end if
    8  continue
       imo=12
    9  continue
#ifndef NOPRINT
#ifdef MP
       if(mype.eq.master) then
#endif
       print 66,fhour,numdyz,jday,nmdtot,ndayr
#ifdef MP
       endif
#endif
   66  format(' sbuvo3 climo hr=',f10.1,
     1        ' numdyz,jday,nmdtot,ndayr=',4i8)
#endif
      end if
      coef=1.655e-6
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      jerr = 0
#ifdef RSM
      open(unit=ko3,file='r_gcmo3',form='formatted',status='old')
#else
      open(unit=ko3,file='gcmo3',form='formatted',status='old')
#endif
      rewind ko3
      do 10 l=1,loz
        read (ko3,15,err=998,end=999) pstr(l)
   15   format(f10.3)
   10 continue
   16 format(1h ,' o3 pressures=',8f10.3)
      do 40 nm=1,12
       do 30 j=1,jmr
         read (ko3,19,err=998,end=999) imon(nm),ilat(j,nm),
     1                                 (o3r(j,l,nm), l=1,10)
         read (ko3,20,err=998,end=999) (o3r(j,l,nm), l=11,loz)
   30  continue
       do 35 j=1,jmr
        do 35 l=1,loz
         o3r(j,l,nm) = o3r(j,l,nm) * coef
   35  continue
   40 continue
#ifdef MP
      endif
#ifndef RSM
      call mpbcasti(jerr,1)
#else
      call rmpbcasti(jerr,1)
#endif
#endif
      if( jerr.ne.0 ) return
#ifdef MP
#ifdef RSM
       call rmpbcastr(pstr,loz)
       call rmpbcastr(o3r,loz*jmr*12)
       call rmpbcasti(imon,12)
       call rmpbcasti(ilat,jmr*12)
#else
      call mpbcastr(pstr,loz)
      call mpbcastr(o3r,loz*jmr*12)
      call mpbcasti(imon,12)
      call mpbcasti(ilat,jmr*12)
#endif
#endif
c      do a linear interpolation in time, where we assume that
c       the ozone data is valid for mid-month
c      monl is the preceeding month, monc for current mo, and
c      monr is the future month..
      monl=imo-1
      monc=imo
      monr=imo+1
      if (monl.lt.1) monl=12
      if (monr.gt.12) monr=1
c
c  difl=number of days beteen mid-months of the current and
c            preceeding mo, difr=same for current and future mo..
c  dell=number of days between current day and mon,
c       delr=same for current day and next month.
c       sign convention as if we were using day of year calculations.
c
      midl=days(monl)/2
      midc=days(monc)/2
      midr=days(monr)/2
      difl=-(days(monl)-midl+midc)
      difr= (days(monc)-midc+midr)
      delday=ida-midc
      if (ida.gt.midc) then
       do 60 j=1,jmr
        do 60 l=1,loz
          o3tmp(j,l)=o3r(j,l,monc)+(o3r(j,l,monr)-o3r(j,l,monc))
     1                * delday/difr
   60  continue
      else if (ida.lt.midc) then
       do 65 j=1,jmr
        do 65 l=1,loz
          o3tmp(j,l)=o3r(j,l,monc)+(o3r(j,l,monl)-o3r(j,l,monc))
     1                * delday/difl
   65  continue
      else if (ida.eq.midc) then
       do 70 j=1,jmr
        do 70 l=1,loz
          o3tmp(j,l)=o3r(j,l,monc)
   70  continue
      end if
#ifndef NOPRINT
#ifdef MP
      if(mype.eq.master) then
#endif
      print 200,imo,ida
#ifdef MP
      endif
#endif
#endif
c
c  flip y-direction for ntrac = 2 for time being
c
      do j=1,jmr
        do l=1,loz
           o3tmp2(j,l)=o3tmp(jmr-j+1,l)
        enddo
      enddo
c
c  linearly interpolate to 5 deg zonal means
c
      jmr1=jmr-1
      do 80 j=1,jmr1
       j1=j*2
       j2=j1+1
       do 80 l=1,loz
        o3out(j1,l)=o3tmp2(j,l)
        o3out(j2,l)=0.5*(o3tmp2(j,l)+o3tmp2(j+1,l))
   80 continue
      do 85 l=1,loz
       o3out(1,l)=o3tmp2(1,l)
       o3out(jmout-1,l)=o3tmp2(jmr,l)
       o3out(jmout,l)=o3tmp2(jmr,l)
   85 continue
   19 format(i2,i4,10f6.2)
   20 format(6x,10f6.2)
  200 format(1h1,'  ozone climatology for month,day=',2i4)
      return
#ifdef MP
#ifdef RSM
#define FAIL rmpabort
#else
#define FAIL mpabort
#endif
#else
#define FAIL abort
#endif
  998 print 988,ko3
      call FAIL
  999 print 989,ko3
      call FAIL
  988 format(1h ,'error reading nasa ozone, unit=',i4)
  989 format(1h ,'e.o.f reading nasa ozone, unit=',i4)
#undef FAIL

      end
