#define JCAP jcap
#define LNT2 lnt2
#define PI pi
      subroutine prtime(id,munth,iyear,ihr,xmin,jd,fjd,
     1                  dlt,alf,r1,slag,solc)
#include <paramodel.h>
#include <comio.h>
      parameter (degrad=180.e0/pi_,hpi=0.5e0*pi_)
      data       sign/1h-/,      sigb/1h /
      data zero,six,sixty,q22855/0.0,6.0,60.0,228.55735/
      save sign,zero,six,sixty,q22855
      dltd=degrad*dlt
      ltd=dltd
      dltm=sixty*(abs(dltd)-abs(float(ltd)))
      ltm=dltm
      dlts=sixty*(dltm-float(ltm))
      dsig=sigb
      if((dltd.lt.zero).and.(ltd.eq.0)) dsig=sign
      halp=six*alf/hpi
      ihalp=halp
      ymin=abs(halp-float(ihalp))*sixty
      iyy=ymin
      asec=(ymin-float(iyy))*sixty
      eqt=q22855*slag
      eqsec=sixty*eqt
#ifndef RSM
#ifndef NOPRINT
      if(iope) print 1004,
     1        id,munth,iyear,ihr,xmin,jd,fjd,r1,halp,ihalp,
     1       iyy,asec,dltd,dsig,ltd,ltm,dlts,eqt,eqsec,slag,solc
 1004 format('0 forecast date',9x,i3,a5,i6,' at',i3,' hrs',f6.2,' mins'/
     1       '  julian day',12x,i8,2x,'plus',f11.6/
     2       '  radius vector',9x,f10.7/
     3       '  right ascension of sun',f12.7,' hrs, or',i4,' hrs',i4,
     4                                 ' mins',f6.1,' secs'/
     5       '  declination of the sun',f12.7,' degs, or',a2,i3,
     6                                 ' degs',i4,' mins',f6.1,' secs'/
     7       '  equation of time',6x,f12.7,' mins, or',f10.2,' secs, or'
     8                           ,f9.6,' radians'/
     9       '  solar constant',8x,f12.7)
#endif
#endif
      return
      end
