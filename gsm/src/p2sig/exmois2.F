#define PI pi
#define PSAT psat
#define RD rd
#define RV rv
#define T0C t0c
      subroutine exmois2(slin,pstar,qin,t,qou,ijdim,kdimi,kdimqi,xlat)
      implicit none
#include <paramodel.h>
c
c
c  this code extrapolates moisture up into stratosphere
c  using exponential decrease up to a specified q at 80 mb,
c  then use lims profile above(see middle atmos program handbook,
c  vol 22,sep 86,ed. j.m.russell ---> j alpert has copy)....
c
c  input - slin  - sigma level 
c        - pstar - surface pressure (cb)
c        - qin   - specific humidity
c        - t     - virtual temperature
c        - xlat  - latitude in radians
c        - ijdim - number of grid points
c        - kdimqi- number of input moisture levels
c        - kdimi - number of output sigma levels
c
c  output - qou - moisture (use lims data abuv 8 cb)
c         - t   - thermodynamic temperature back into t
c
      integer ijdim
      integer kdimi,kdimqi
      real slin(kdimi)
      real pstar(ijdim)
      real qin(ijdim,kdimqi)
      real t(ijdim,kdimi)
      real xlat(ijdim)
      real qou(ijdim,kdimi)
c
      real qout(kdimi)
      real qsat(kdimi),prs(kdimi)
c
      integer i,j,k,ij
      integer jdx,ldry,kdy
      real pmin,pmm,xlnpm
      real expon,es,dlat,rh,dx,q8,q1,xlnqm,dy
      real qqout,xlnpq,xlnqq,xlamb
c
c  qq stores lims (satellite-1986-see j.alpert) observed h2o above
c  80 mb...(k,j)-k=1,8 for 80,70,..,10 mb
c  j=1,19 for every 10 deg lat (starting at n.pole)
c
      real qq(8,19)
      save qq
      data qq /32*5.  ,  4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7,
     1 4.,3.85,3.9,4.,4.25,4.5,4.5,4.6, 3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,
     2           2.1,2.,2.,2.5,3.,3.9,4.,4. ,
     3 2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.,2.,2.,2.3,2.7,3.25,3.45,3.75,
     4 2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.1,2.,2.,2.5,3.,3.9,4.,4. ,
     5           3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,
     6      4.,3.85,3.9,4.,4.25,4.5,4.5,4.6 ,
     7           4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7 , 32*5.  /
c
c  first convert mdl virtual temp to thermodynmic temp
c
      pmin = 8.0
      pmm  = 1.0
      xlnpm =  log(pmin)
c
c  get correct units of lims moisture
c
      do j=1,19
        do k=1,8
          qq(k,j) = qq(k,j) * 1.e-6
        enddo
      enddo
c
      do i=1,ijdim
c
c  get latitude in degrees..
c
      dlat = xlat(i) * 180.0/pi_
c
c  model temp is virtual (since hydrostatic from hgts initially),
c
      do k=1,kdimqi
        if (qin(i,k).le.0.0) go to 2
        t(i,k)=t(i,k)*(1.0+qin(i,k))/(1.0+rv_/rd_*qin(i,k))
   2    continue
      enddo
c
c  get lyr pressure and then saturated moisture
c
      do k=1,kdimi
        prs(k) = pstar(i) * slin(k)
      enddo
c
c  compute saturation specific humidity(dimensionless) from
c  temperature t (deg k) and pressure (cb)
c  conversion to specific humidity follows from definition
c
      do k=1,kdimi
        expon = 7.50*(t(i,k)-t0c_)/((t(i,k)-t0c_)+237.30)
        es = psat_*1.e-2 * 10.0**expon
        qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)
      enddo
c
c  limit moisture in lowest levh layers--rh le 1. but ge .15
c          the latter avoids the negative q problem...
c
      do k=1,kdimqi
        qout(k) = qin(i,k)
        rh = qin(i,k) / qsat(k)
      enddo
c
c  obtain 8.cb (k=1) and 1.cb (k=8) valu by horiz interpolation
c
      jdx = (90.-dlat)/10.+1.
      dx = (90.-dlat)/10.+1.-jdx
      q8  = qq(1,jdx)*(1.-dx)+qq(1,jdx+1)*dx
      q1  = qq(8,jdx)*(1.-dx)+qq(8,jdx+1)*dx
      ldry = kdimqi + 1
c 
      if( ldry .le. kdimi ) then
c
        xlnqm =  log(q8)
c
c  extrapolate moisture to min valu(q8) at pressure pmin(80 mb)
c   use exponential decrease from layer levh----
c   i.e.  q=q(levh)*(p/p(levh)) ** xlamb
c   where xlamb is computed to fit q8 at pmin and q,p at levh.
c
        qqout = qin(i,kdimqi)
        rh = qqout / qsat(kdimqi)
        if (rh.le..15) qqout = .15 * qsat(kdimqi)
        if (rh.gt. 1.0) qqout = qsat(kdimqi)
        xlnpq =  log(prs(kdimqi))
        xlnqq =  log(qqout)
        xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq)
        do k=ldry,kdimi
          if (prs(k).ge.pmin) go to 21
          if (prs(k).ge.pmm ) go to 22
c
c  above pmm(10 mb) use constant value from table
c
          qout(k) = q1
          go to 13
   21     continue
          qout(k) = qout(kdimqi)*(prs(k)/prs(kdimqi)) ** xlamb
          go to 13
   22     continue
c
c  above 8 cb so complete linear interp from table..
c
          kdy = 9. - prs(k)
          dy = 9. - prs(k) - kdy
          qout(k) = qq(kdy,jdx  )*(1.-dy)*(1.-dx)+qq(kdy+1,jdx+1)*dy*dx
     1             +qq(kdy,jdx+1)*(1.-dy)*dx+qq(kdy+1,jdx)*dy*(1.-dx)
   13     rh = qout(k) / qsat(k)
          if(rh.gt.1.) qout(k)=qsat(k)
        enddo
      endif
c
c  store extrapolated moisture
c
      do k=1,kdimi
        qou(i,k) = qout(k)
      enddo
c
c  temp. is restored to virtual temp.
c
      do k=1,kdimqi
        if (qin(i,k).gt.0.) then
          t(i,k) = t(i,k) * (1. + rv_/rd_ * qin(i,k))
     1                   /(1. + qin(i,k))
        endif
      enddo
c
      enddo
c
      return
      end
