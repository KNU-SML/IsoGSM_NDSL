#define JCAP jcap
#define LNT2 lnt2
#define PI pi
      subroutine fcstim(fhour,imon,iday,iztim,jdnmc,fjdnmc,rlag,year,
     1                  rsin1,rcos1,rcos2,jd,fjd)
#include <paramodel.h>
#include <comio.h>
      parameter (tpi=2.e0*pi_)
      dimension jmon(12)
      data jmon/31,28,31,30,31,30,31,31,30,31,30,31/
      data two/2.e0/
      save jmon,two
c....    first get number of days since beginning of year (no leap yrs)
      nnday =0
      imo = imon - 1
      if (imo.gt.0) then
        do 150 i=1,imo
          nnday = nnday + jmon(i)
  150   continue
      endif
      nnday = nnday + iday
#ifndef RSM
#ifndef NOPRINT
      if(iope) print 1002,nnday
 1002 format(1h ,'*************** nnday of year = ',i4,'******')
#endif
#endif
c....    get number of days into fcst (dyfcst)
c....   following two cards changed on 10 apr 86 to fix slight error
c         in solar declination calc if initial hr not 00z or 12z.....
      dayini = nnday + float(iztim)/24.e0
      soltim = fhour + iztim
c...     reset to 24 hour clock
      fday = soltim / 24.e0
c>yh  soltim = soltim - int(fday) * 24.e0
      dyfcst = nnday + fday
#ifndef RSM
#ifndef NOPRINT
      if(iope) print 1003,jdnmc,fjdnmc,fhour,dayini,dyfcst
 1003 format(1h0,'from heatl3 jdnmc etc',i9,2x,4(2x,f6.2))
#endif
#endif
      rang=tpi*(dyfcst-rlag)/year
      rsin1=sin(rang)
      rcos1=cos(rang)
      rcos2=cos(two*rang)
c....    update the julian date (initial in jdnmc,fjdnmc)
      dyinc = dyfcst - dayini
      idyin = dyinc
      fdyin = dyinc - idyin
      jd = jdnmc + idyin
      fjd = fjdnmc + fdyin
c.......need to reset if fraction (fjd) gt 1.
      ifjd = fjd
      if (ifjd.gt.0) then
        jd = jd + ifjd
        fjd = fjd - ifjd
      endif
      return
      end
