#define CP cp
#define JCAP jcap
#define LEVS levs
#define LNT2 lnt2
#define RD rd
      subroutine gcljms (si)
#include <paramodel.h>
      parameter (kdim=levs_, kdimp=kdim+1)
#include <comcd1.h>
#include <comio.h>
      dimension si(kdimp), ppptop(4,2)
c --- pressure limits for sfc and top of each cloud domain (l,m,h)
c     in mb, model layers for cld tops are l=7,m=11,h=15 at low
c     latitudes and l= ,m= ,h=  , at pole region.
c....     ptop above h changed from 150 to 100, cause
cc    data ppptop /1050.,642.,350.,150., 1050.,750.,500.,150./
c          code was truncating tops of convective clouds
      data ppptop /1050.,642.,350.,100., 1050.,750.,500.,100./
c
      rocp = rd_ / cp_
c --- inverson type cld critical value-istrat=0
cyh94 clapse = -0.055e0
      clapse = -0.06e0
c --- inverson type cld critical value-istrat=1
      clapkc = -0.05e0
c....critical dtheta/dp for ocean stratus(wgt varies 0 to 1
c                linearly from clapse to clpse)
      dclps = -0.01e0
      clpse = clapkc + dclps
      cvtop = 400.0e0
      pstrt = 800.0e0
c --- low cld bottom (at sigma=0.95) and top sigma level
      do 5 k=1,kdim
        kk=k
        if (si(kk) .le. 0.95e0) go to 10
   5  continue
  10  klowb = kk - 1
      silow = ppptop(2,1) * 1.0e-3
      do 20 k=1,kdim
        kk=k
        if (si(kk) .lt. silow) go to 30
  20  continue
  30  klowt = kk
c --- presure limit at sfc and at top of cloud domains (l,m,h) in mb
      do 40 j=1,2
      do 40 i=1,4
       ptopc(i,j) = ppptop(i,j)
  40  continue
c --- l cld vertical vel adj boundaries
c
#ifndef CLD1
cc      vvcld(1) =  0.0003e0
cc      vvcld(2) = -0.0005e0
c  turned off vertical motion check
      vvcld(1) =  100.
      vvcld(2) =  100.
c
c  changed by mk 5/5/98
c
c     vvcld(1) =  0.0006
c     vvcld(2) =  0.0006
#else
c  turned off vertical motion check
      vvcld(1) =  100.
      vvcld(2) =  100.
#endif
c
      crhrh = 0.60e0
c--- compute llyr--which is topmost non cld(low) layer, for stratiform
      xthk = 0.e0
c....   default llyr
      kl = kdimp
c....   topmost noncloud layer will be the one at or above lowest
c         0.1 of the atmosphere..
      do 202 k=1,kdim
c       xthk = xthk + si(k) - si(k+1)
c       if (xthk.lt.0.1e0) go to 202
        kl = k
c       go to 204
        if (si(k).lt.0.9e0) go to 204
  202 continue
  204 llyr = kl-1
#ifndef NOPRINT
      if(iope) print 205,llyr,klowb
  205 format(1h ,'-------llyr,klowb =',2i5)
#endif
      return
      end
