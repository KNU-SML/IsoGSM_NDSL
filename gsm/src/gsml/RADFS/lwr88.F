#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
#include <define.h>
      subroutine lwr88(ipts,heatra,grnflx,topflx,
#ifdef CLR
     1                 grnfx0,topfx0,
#endif
     1                 press,temp,rh2o,qo3,cldfac,
     2                 camt,nclds,ktop,kbtm)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c-paramodel
c     subroutine lwr88 computes temperature-corrected co2 transmission
c   functions and also computes the pressure grid and layer optical
c   paths.
c          inputs:                (common blocks)
c      cldfac                          cldcom
c      press,temp,rh2o,qo3             radisw
c      camt,nclds,ktop,kbtm            radisw
c      co251,co258,cdt51,cdt58         co2bd3
c      c2d51,c2d58,co2m51,co2m58       co2bd3
c      cdtm51,cdtm58,c2dm51,c2dm58     co2bd3
c      stemp,gtemp                     co2bd3
c      co231,co238,cdt31,cdt38         co2bd2
c      c2d31,c2d38                     co2bd2
c      co271,co278,cdt71,cdt78         co2bd4
c      c2d71,c2d78                     co2bd4
c      betinw                          bdwide
c          outputs:
c      heatra,grnflx,topflx            lwout
c          called by:
c      radmn or input routine of model
c          calls:
c      fst88
c
#include <hcon.h>
#include <rdparm.h>
#include <co2dta.h>
#include <rnddta.h>
c
      dimension press(imbx,lp1),temp(imbx,lp1),rh2o(imbx,l),qo3(imbx,l)
      dimension cldfac(imbx,lp1,lp1),camt(imbx,lp1)
      dimension nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)
      dimension heatra(imbx,l),grnflx(imax),topflx(imax)
#ifdef CLR
      dimension grnfx0(imax),topfx0(imax)
#endif
      dimension delp2(imbx,l)
c
      dimension qh2o(imbx,l),t(imbx,lp1)
      dimension p(imbx,lp1),delp(imbx,l)
      dimension co21(imbx,lp1,lp1),co2nbl(imbx,l)
      dimension co2sp1(imbx,lp1),co2sp2(imbx,lp1)
      dimension var1(imbx,l),var2(imbx,l),var3(imbx,l),var4(imbx,l)
      dimension cntval(imbx,lp1)
      dimension toto3(imbx,lp1),tphio3(imbx,lp1),totphi(imbx,lp1)
      dimension totvo2(imbx,lp1),emx1(imax),emx2(imax),empl(imbx,llp1)
c
      dimension co2r(imbx,lp1),dift(imbx,lp1)
      dimension co2r1(imbx,lp1),dco2d1(imbx,lp1)
      dimension d2cd21(imbx,lp1),d2cd22(imbx,lp1)
      dimension co2r2(imbx,lp1),dco2d2(imbx,lp1)
      dimension co2mr(imbx,l),co2md(imbx,l),co2m2d(imbx,l)
      dimension tdav(imbx,lp1),tstdav(imbx,lp1),
     1 vv(imbx,l),vsum3(imbx,lp1),vsum1(imax),vsum2(imax)
      dimension a1(imax),a2(imax)
      dimension dco2dt(imbx,lp1),d2cdt2(imbx,lp1)
c
      dimension texpsl(imbx,lp1),tlsqu(imbx,lp1)
      dimension vsum4(imbx,l)
      equivalence (vsum3,tlsqu,texpsl)
      equivalence (vv,vsum4)
c
c****compute flux pressures (p) and differences (delp2,delp)
c****compute flux level temperatures (t) and continuum temperature
c    corrections (texpsl)
      do 103 k=2,l
      do 103 i=1,ipts
      p(i,k)=haf*(press(i,k-1)+press(i,k))
      t(i,k)=haf*(temp(i,k-1)+temp(i,k))
103   continue
      do 105 i=1,ipts
      p(i,1)=zero
      p(i,lp1)=press(i,lp1)
      t(i,1)=temp(i,1)
      t(i,lp1)=temp(i,lp1)
105   continue
      do 107 k=1,l
      do 107 i=1,ipts
      delp2(i,k)=p(i,k+1)-p(i,k)
      delp(i,k)=one/delp2(i,k)
107   continue
c****compute argument for cont.temp.coeff.
c    (this is 1800.(1./temp-1./296.))
cmk   if( temp(ipts,lp1).lt.1. ) close(90)
      do 125 k=1,lp1
      do 125 i=1,ipts
      texpsl(i,k)=h18e3/temp(i,k)-h6p08108
c...then take exponential
      texpsl(i,k)=exp(texpsl(i,k))
125   continue
c***compute optical paths for h2o and o3, using the diffusivity
c   approximation for the angular integration (1.66). obtain the
c   unweighted values(var1,var3) and the weighted values(var2,var4).
c   the quantities h3m4(.0003) and h3m3(.003) appearing in the var2 and
c   var4 expressions are the approximate voigt corrections for h2o and
c   o3,respectively.
c
      do 131 k=1,l
      do 131 i=1,ipts
      qh2o(i,k)=rh2o(i,k)*diffctr
c---vv is the layer-mean pressure (in atm),which is not the same as
c   the level pressure (press)
      vv(i,k)=haf*(p(i,k+1)+p(i,k))*p0inv
      var1(i,k)=delp2(i,k)*qh2o(i,k)*ginv
      var3(i,k)=delp2(i,k)*qo3(i,k)*diffctr*ginv
      var2(i,k)=var1(i,k)*(vv(i,k)+h3m4)
      var4(i,k)=var3(i,k)*(vv(i,k)+h3m3)
c  compute optical path for the h2o continuum, using roberts coeffs.
c  (betinw),and temp. correction (texpsl). the diffusivity factor
c  (which cancels out in this expression) is assumed to be 1.66. the
c  use of the diffusivity factor has been shown to be a significant
c  source of error in the continuum calcs.,but the time penalty of
c  an angular integration is severe.
c
      cntval(i,k)=texpsl(i,k)*rh2o(i,k)*var2(i,k)*betinw/
     1             (rh2o(i,k)+rath2omw)
131   continue
c   compute summed optical paths for h2o,o3 and continuum
      do 201 i=1,ipts
      totphi(i,1)=zero
      toto3(i,1)=zero
      tphio3(i,1)=zero
      totvo2(i,1)=zero
201   continue
      do 203 k=2,lp1
      do 203 i=1,ipts
      totphi(i,k)=totphi(i,k-1)+var2(i,k-1)
      toto3(i,k)=toto3(i,k-1)+var3(i,k-1)
      tphio3(i,k)=tphio3(i,k-1)+var4(i,k-1)
      totvo2(i,k)=totvo2(i,k-1)+cntval(i,k-1)
203   continue
c---emx1 is the additional pressure-scaled mass from press(l) to
c   p(l). it is used in nearby layer and emiss calculations.
c---emx2 is the additional pressure-scaled mass from press(l) to
c   p(lp1). it is used in calculations between flux levels l and lp1.
c
      do 801 i=1,ipts
      emx1(i)=qh2o(i,l)*press(i,l)*(press(i,l)-p(i,l))*gp0inv
      emx2(i)=qh2o(i,l)*press(i,l)*(p(i,lp1)-press(i,l))*gp0inv
801   continue
c---empl is the pressure scaled mass from p(k) to press(k) (index 2-lp1)
c   or to press(k+1) (index lp2-ll)
      do 811 k=1,l
      do 811 i=1,ipts
      empl(i,k+1)=qh2o(i,k)*p(i,k+1)*(p(i,k+1)-press(i,k))*gp0inv
811   continue
      do 812 k=1,lm1
      do 812 i=1,ipts
      empl(i,k+lp1)=qh2o(i,k+1)*p(i,k+1)*(press(i,k+1)-p(i,k+1))*gp0inv
812   continue
      do 821 i=1,ipts
      empl(i,1)=var2(i,l)
      empl(i,llp1)=empl(i,ll)
821   continue
c***compute weighted temperature (tdav) and pressure (tstdav) integrals
c   for use in obtaining temp. difference bet. sounding and std.
c   temp. sounding (dift)
      do 161 i=1,ipts
      tstdav(i,1)=zero
      tdav(i,1)=zero
161   continue
      do 162 k=1,lp1
      do 162 i=1,ipts
      vsum3(i,k)=temp(i,k)-stemp(k)
162   continue
      do 163 k=1,l
      do 165 i=1,ipts
      vsum2(i)=gtemp(k)*delp2(i,k)
      vsum1(i)=vsum2(i)*vsum3(i,k)
      tstdav(i,k+1)=tstdav(i,k)+vsum2(i)
      tdav(i,k+1)=tdav(i,k)+vsum1(i)
165   continue
163   continue
c
c****evaluate coefficients for co2 pressure interpolation (a1,a2)
      do 171 i=1,ipts
      a1(i)=(press(i,lp1)-p0xzp8)/p0xzp2
      a2(i)=(p0-press(i,lp1))/p0xzp2
171   continue
c***perform co2 pressure interpolation on all inputted transmission
c   functions and temp. derivatives
c---successively computing co2r,dco2dt and d2cdt2 is done to save
c   storage (at a slight loss in computation time)
      do 184 k=1,lp1
      do 184 i=1,ipts
        co2r1(i,k)=a1(i)*co231(k)+a2(i)*co238(k)
        d2cd21(i,k)=h1m3*(a1(i)*c2d31(k)+a2(i)*c2d38(k))
        dco2d1(i,k)=h1m2*(a1(i)*cdt31(k)+a2(i)*cdt38(k))
        co2r2(i,k)=a1(i)*co271(k)+a2(i)*co278(k)
        d2cd22(i,k)=h1m3*(a1(i)*c2d71(k)+a2(i)*c2d78(k))
        dco2d2(i,k)=h1m2*(a1(i)*cdt71(k)+a2(i)*cdt78(k))
184   continue
      do 190 k=1,l
      do 190 i=1,ipts
        co2mr(i,k)=a1(i)*co2m51(k)+a2(i)*co2m58(k)
        co2md(i,k)=h1m2*(a1(i)*cdtm51(k)+a2(i)*cdtm58(k))
        co2m2d(i,k)=h1m3*(a1(i)*c2dm51(k)+a2(i)*c2dm58(k))
190   continue
c***compute co2 temperature interpolations for all bands,using dift
c
c   the case where k=1 is handled first. we are now replacing
c   3-dimensional arrays by 2-d arrays, to save space. thus this
c   calculation is for (i,kp,1)
      do 211 kp=2,lp1
      do 211 i=1,ipts
      dift(i,kp)=tdav(i,kp)/tstdav(i,kp)
211   continue
      do 212 i=1,ipts
      co21(i,1,1)=1.0
      co2sp1(i,1)=1.0
      co2sp2(i,1)=1.0
212   continue
      do 215 kp=2,lp1
      do 215 i=1,ipts
c---calculations for kp>1 for k=1
      co2r(i,kp)=a1(i)*co251(kp,1)+a2(i)*co258(kp,1)
      dco2dt(i,kp)=h1m2*(a1(i)*cdt51(kp,1)+a2(i)*cdt58(kp,1))
      d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(kp,1)+a2(i)*c2d58(kp,1))
      co21(i,kp,1)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+
     1             haf*dift(i,kp)*d2cdt2(i,kp))
c---calculations for (effectively) kp=1,k>kp. these use the
c   same value of dift due to symmetry
      co2r(i,kp)=a1(i)*co251(1,kp)+a2(i)*co258(1,kp)
      dco2dt(i,kp)=h1m2*(a1(i)*cdt51(1,kp)+a2(i)*cdt58(1,kp))
      d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(1,kp)+a2(i)*c2d58(1,kp))
      co21(i,1,kp)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+
     1             haf*dift(i,kp)*d2cdt2(i,kp))
215   continue
c   the transmission functions used in spa88 may be computed now.
c---(in the 250 loop,dift really should be (i,1,k), but dift is
c    invariant with respect to k,kp,and so (i,1,k)=(i,k,1))
      do 250 k=2,lp1
      do 250 i=1,ipts
      co2sp1(i,k)=co2r1(i,k)+dift(i,k)*(dco2d1(i,k)+haf*dift(i,k)*
     1 d2cd21(i,k))
      co2sp2(i,k)=co2r2(i,k)+dift(i,k)*(dco2d2(i,k)+haf*dift(i,k)*
     1 d2cd22(i,k))
250   continue
c
c   next the case when k=2...l
      do 220 k=2,l
      do 222 kp=k+1,lp1
      do 222 i=1,ipts
      dift(i,kp)=(tdav(i,kp)-tdav(i,k))/
     1              (tstdav(i,kp)-tstdav(i,k))
      co2r(i,kp)=a1(i)*co251(kp,k)+a2(i)*co258(kp,k)
      dco2dt(i,kp)=h1m2*(a1(i)*cdt51(kp,k)+a2(i)*cdt58(kp,k))
      d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(kp,k)+a2(i)*c2d58(kp,k))
      co21(i,kp,k)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+
     1             haf*dift(i,kp)*d2cdt2(i,kp))
      co2r(i,kp)=a1(i)*co251(k,kp)+a2(i)*co258(k,kp)
      dco2dt(i,kp)=h1m2*(a1(i)*cdt51(k,kp)+a2(i)*cdt58(k,kp))
      d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(k,kp)+a2(i)*c2d58(k,kp))
      co21(i,k,kp)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+
     1             haf*dift(i,kp)*d2cdt2(i,kp))
222   continue
220   continue
c   finally the case when k=kp,k=2..lp1
      do 206 k=2,lp1
      do 206 i=1,ipts
      dift(i,k)=haf*(vsum3(i,k)+vsum3(i,k-1))
      co2r(i,k)=a1(i)*co251(k,k)+a2(i)*co258(k,k)
      dco2dt(i,k)=h1m2*(a1(i)*cdt51(k,k)+a2(i)*cdt58(k,k))
      d2cdt2(i,k)=h1m3*(a1(i)*c2d51(k,k)+a2(i)*c2d58(k,k))
      co21(i,k,k)=co2r(i,k)+dift(i,k)*(dco2dt(i,k)+
     1             haf*dift(i,k)*d2cdt2(i,k))
206   continue
c--- we arent doing nbl tfs on the 100 cm-1 bands .
      do 260 k=1,l
      do 260 i=1,ipts
      co2nbl(i,k)=co2mr(i,k)+vsum3(i,k)*(co2md(i,k)+haf*
     1 vsum3(i,k)*co2m2d(i,k))
260   continue
c***compute temp. coefficient based on t(k) (see ref.2)
      do 264 k=1,lp1
      do 264 i=1,ipts
      if (t(i,k).le.h25e2) then
         tlsqu(i,k)=b0+(t(i,k)-h25e2)*
     1                      (b1+(t(i,k)-h25e2)*
     2                   (b2+b3*(t(i,k)-h25e2)))
      else
         tlsqu(i,k)=b0
      endif
264   continue
c***apply to all co2 tfs
      do 280 k=1,lp1
      do 282 kp=1,lp1
      do 282 i=1,ipts
      co21(i,kp,k)=co21(i,kp,k)*(one-tlsqu(i,kp))+tlsqu(i,kp)
282   continue
280   continue
      do 284 k=1,lp1
      do 286 i=1,ipts
      co2sp1(i,k)=co2sp1(i,k)*(one-tlsqu(i,1))+tlsqu(i,1)
      co2sp2(i,k)=co2sp2(i,k)*(one-tlsqu(i,1))+tlsqu(i,1)
286   continue
284   continue
      do 288 k=1,l
      do 290 i=1,ipts
      co2nbl(i,k)=co2nbl(i,k)*(one-tlsqu(i,k))+tlsqu(i,k)
290   continue
288   continue
      call fst88(ipts,heatra,grnflx,topflx,
#ifdef CLR
     1           grnfx0,topfx0,
#endif
     1           qh2o,press,p,delp,delp2,temp,t,
     2           cldfac,nclds,ktop,kbtm,camt,
     3           co21,co2nbl,co2sp1,co2sp2,
     4           var1,var2,var3,var4,cntval,
     5           toto3,tphio3,totphi,totvo2,
     6           emx1,emx2,empl)
      return
      end
