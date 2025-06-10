#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
#include <define.h>
      subroutine fst88(ipts,heatra,grnflx,topflx,
#ifdef CLR
     1                 grnfx0,topfx0,
#endif
     1                 qh2o,press,p,delp,delp2,temp,t,
     2                 cldfac,nclds,ktop,kbtm,camt,
     3                 co21,co2nbl,co2sp1,co2sp2,
     4                 var1,var2,var3,var4,cntval,
     5                 toto3,tphio3,totphi,totvo2,
     6                 emx1,emx2,empl)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c-paramodel
c
#include <hcon.h>
#include <rdparm.h>
#include <rnddta.h>
#include <tabcom.h>
c
      dimension qh2o(imbx,lp1),press(imbx,lp1)
      dimension p(imbx,lp1),delp(imbx,l),delp2(imbx,l),temp(imbx,lp1)
      dimension t(imbx,lp1),cldfac(imbx,lp1,lp1),camt(imbx,lp1)
      dimension nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)
      dimension co21(imbx,lp1,lp1),co2nbl(imbx,l)
      dimension co2sp1(imbx,lp1),co2sp2(imbx,lp1)
      dimension var1(imbx,l),var2(imbx,l),var3(imbx,l),var4(imbx,l)
      dimension cntval(imbx,lp1)
      dimension heatra(imbx,l),grnflx(imax),topflx(imax)
c
#ifdef CLR
      dimension heatr0(imbx,l),flxnt0(imbx,lp1)
      dimension grnfx0(imax),topfx0(imax),gxcts0(imax),flx1e0(imax)
      dimension excts0(imbx,l),ctso30(imbx,l),cts0(imbx,l)
      dimension flx0(imbx,lp1)
#endif
c
      dimension gxcts(imax),flx1e1(imax)
      dimension avephi(imbx,lp1),emiss(imbx,lp1),emissb(imbx,lp1)
c
      dimension toto3(imbx,lp1),tphio3(imbx,lp1),totphi(imbx,lp1)
      dimension totvo2(imbx,lp1),emx1(imax),emx2(imax),empl(imbx,llp1)
c
      dimension excts(imbx,l),ctso3(imbx,l),cts(imbx,l),e1flx(imbx,lp1)
      dimension co2sp(imbx,lp1),to3spc(imbx,l),to3sp(imbx,lp1)
      dimension oss(imbx,lp1),css(imbx,lp1),ss1(imbx,lp1),ss2(imbx,lp1),
     1 tc(imbx,lp1),dtc(imbx,lp1)
      dimension sorc(imbx,lp1,nbly),csour(imbx,lp1)
ccc
      dimension        avvo2(imbx,lp1),heatem(imbx,lp1),
     1                 over1d(imbx,lp1),
     1                 to31d(imbx,lp1),cont1d(imbx,lp1),
     2                 avmo3(imbx,lp1),avpho3(imbx,lp1),
     2                 c(imbx,llp1),c2(imbx,llp1)
      dimension       itop(imax),ibot(imax),indtc(imax)
      dimension
     4 delptc(imax),ptop(imax),pbot(imax),ftop(imax),
     5 fbot(imax) ,emspec(imbx,2)
c---dimension of variables equivalenced to those in vtemp---
      dimension vtmp3(imbx,lp1),dsorc(imbx,lp1)
      dimension alp(imbx,llp1),csub(imbx,llp1),csub2(imbx,llp1)
      dimension fac1(imbx,lp1)
      dimension delpr1(imbx,lp1),delpr2(imbx,lp1)
      dimension emisdg(imbx,lp1),contdg(imbx,lp1),to3dg(imbx,lp1)
      dimension flxnet(imbx,lp1)
      dimension ixo(imbx,lp1)
      dimension vsum1(imbx,lp1)
      dimension flxthk(imbx,lp1)
      dimension z1(imbx,lp1)
c---dimension of variables passed to other subroutines---
c   (and not found in common blocks)
      dimension e1cts1(imbx,lp1),e1cts2(imbx,l)
      dimension e1ctw1(imbx,lp1),e1ctw2(imbx,l)
      dimension emd(imbx,llp1),tpl(imbx,llp1)
c   it is possible to equivalence emd,tpl to the above variables,
c   as they get called at different times
      dimension fxo(imbx,lp1),dt(imbx,lp1)
      dimension fxoe2(imbx,lp1),dte2(imbx,lp1)
      dimension fxosp(imbx,2),dtsp(imbx,2)
c
c     dimension of local variables
      dimension rlog(imbx,l),flx(imbx,lp1)
      dimension totevv(imbx,lp1),cnttau(imbx,lp1)
c
      equivalence (alp,c,csub),(csub2,c2)
      equivalence (fac1,dsorc,over1d,delpr2,flxnet)
      equivalence (delpr1,heatem)
      equivalence (ixo,avvo2,flxthk,to3dg)
      equivalence (z1,avmo3,contdg)
      equivalence (emisdg,vsum1,avpho3)
      equivalence (emd(1,1),e1cts1(1,1)),(emd(1,lp2),e1cts2(1,1))
      equivalence (tpl(1,1),e1ctw1(1,1)),(tpl(1,lp2),e1ctw2(1,1))
c
c          first section is table lookup for source function and
c     derivative (b and db/dt).also,the nlte co2 source function
c     is obtained
c
c---in calcs. below, decrementing the index by 9
c   accounts for the tables beginning at t=100k.
c   at t=100k.
      do 101 k=1,lp1
      do 101 i=1,ipts
c---temp. indices for e1,source
      vtmp3(i,k)=aint(temp(i,k)*hp1)
      fxo(i,k)=vtmp3(i,k)-9.
      dt(i,k)=temp(i,k)-ten*vtmp3(i,k)
c---integer index for source (used immediately)
      ixo(i,k)=fxo(i,k)
101   continue
      do 103 k=1,l
      do 103 i=1,ipts
c---temp. indices for e2 (kp=1 layer not used in flux calculations)
      vtmp3(i,k)=aint(t(i,k+1)*hp1)
      fxoe2(i,k)=vtmp3(i,k)-9.
      dte2(i,k)=t(i,k+1)-ten*vtmp3(i,k)
103   continue
c---special case to handle kp=lp1 layer and special e2 calcs.
      do 105 i=1,ipts
      fxoe2(i,lp1)=fxo(i,l)
      dte2(i,lp1)=dt(i,l)
      fxosp(i,1)=fxoe2(i,lm1)
      fxosp(i,2)=fxo(i,lm1)
      dtsp(i,1)=dte2(i,lm1)
      dtsp(i,2)=dt(i,lm1)
105   continue
c
c---source function for combined band 1
      do 4114 i=1,ipts
      do 4114 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),1)
        dsorc(i,k)=dsrce(ixo(i,k),1)
4114   continue
      do 4112 k=1,lp1
      do 4112 i=1,ipts
      sorc(i,k,1)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4112   continue
c---source function for combined band 2
      do 4214 i=1,ipts
      do 4214 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),2)
        dsorc(i,k)=dsrce(ixo(i,k),2)
4214   continue
      do 4212 k=1,lp1
      do 4212 i=1,ipts
      sorc(i,k,2)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4212   continue
c---source function for combined band 3
      do 4314 i=1,ipts
      do 4314 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),3)
        dsorc(i,k)=dsrce(ixo(i,k),3)
4314   continue
      do 4312 k=1,lp1
      do 4312 i=1,ipts
      sorc(i,k,3)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4312   continue
c---source function for combined band 4
      do 4414 i=1,ipts
      do 4414 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),4)
        dsorc(i,k)=dsrce(ixo(i,k),4)
4414   continue
      do 4412 k=1,lp1
      do 4412 i=1,ipts
      sorc(i,k,4)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4412   continue
c---source function for combined band 5
      do 4514 i=1,ipts
      do 4514 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),5)
        dsorc(i,k)=dsrce(ixo(i,k),5)
4514   continue
      do 4512 k=1,lp1
      do 4512 i=1,ipts
      sorc(i,k,5)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4512   continue
c---source function for combined band 6
      do 4614 i=1,ipts
      do 4614 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),6)
        dsorc(i,k)=dsrce(ixo(i,k),6)
4614   continue
      do 4612 k=1,lp1
      do 4612 i=1,ipts
      sorc(i,k,6)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4612   continue
c---source function for combined band 7
      do 4714 i=1,ipts
      do 4714 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),7)
        dsorc(i,k)=dsrce(ixo(i,k),7)
4714   continue
      do 4712 k=1,lp1
      do 4712 i=1,ipts
      sorc(i,k,7)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4712   continue
c---source function for combined band 8
      do 4814 i=1,ipts
      do 4814 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),8)
        dsorc(i,k)=dsrce(ixo(i,k),8)
4814   continue
      do 4812 k=1,lp1
      do 4812 i=1,ipts
      sorc(i,k,8)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4812   continue
c---source function for band 9 (560-670 cm-1)
      do 4914 i=1,ipts
      do 4914 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),9)
        dsorc(i,k)=dsrce(ixo(i,k),9)
4914   continue
      do 4912 k=1,lp1
      do 4912 i=1,ipts
      sorc(i,k,9)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
4912   continue
c---source function for band 10 (670-800 cm-1)
      do 5014 i=1,ipts
      do 5014 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),10)
        dsorc(i,k)=dsrce(ixo(i,k),10)
5014  continue
      do 5012 k=1,lp1
      do 5012 i=1,ipts
      sorc(i,k,10)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
5012   continue
c---source function for band 11 (800-900 cm-1)
      do 5114 i=1,ipts
      do 5114 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),11)
        dsorc(i,k)=dsrce(ixo(i,k),11)
5114   continue
      do 5112 k=1,lp1
      do 5112 i=1,ipts
      sorc(i,k,11)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
5112   continue
c---source function for band 12 (900-990 cm-1)
      do 5214 i=1,ipts
      do 5214 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),12)
        dsorc(i,k)=dsrce(ixo(i,k),12)
5214   continue
      do 5212 k=1,lp1
      do 5212 i=1,ipts
      sorc(i,k,12)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
5212   continue
c---source function for band 13 (990-1070 cm-1)
      do 5314 i=1,ipts
      do 5314 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),13)
        dsorc(i,k)=dsrce(ixo(i,k),13)
5314   continue
      do 5312 k=1,lp1
      do 5312 i=1,ipts
      sorc(i,k,13)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
5312   continue
c---source function for band 14 (1070-1200 cm-1)
      do 5414 i=1,ipts
      do 5414 k=1,lp1
        vtmp3(i,k)=source(ixo(i,k),14)
        dsorc(i,k)=dsrce(ixo(i,k),14)
5414   continue
      do 5412 k=1,lp1
      do 5412 i=1,ipts
      sorc(i,k,14)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
5412   continue
c
c        the following subroutine obtains nlte source function for co2
c
c
c     call nlte
c
c
c---obtain special source functions for the 15 um band (csour)
c   and the window region (ss1)
      do 131 k=1,lp1
      do 131 i=1,ipts
      ss1(i,k)=sorc(i,k,11)+sorc(i,k,12)+sorc(i,k,14)
131   continue
      do 143 k=1,lp1
      do 143 i=1,ipts
      csour(i,k)=sorc(i,k,9)+sorc(i,k,10)
143   continue
c
c---compute temp**4 (tc) and vertical temperature differences
c   (oss,css,ss2,dtc). all these will be used later in flux computa-
c   tions.
c
      do 901 k=1,lp1
      do 901 i=1,ipts
      tc(i,k)=(temp(i,k)*temp(i,k))**2
901   continue
      do 903 k=1,l
      do 903 i=1,ipts
      oss(i,k+1)=sorc(i,k+1,13)-sorc(i,k,13)
      css(i,k+1)=csour(i,k+1)-csour(i,k)
      dtc(i,k+1)=tc(i,k+1)-tc(i,k)
      ss2(i,k+1)=ss1(i,k+1)-ss1(i,k)
903   continue
c
c
c---the followimg is a drastic rewrite of the radiation code to
c    (largely) eliminate three-dimensional arrays. the code works
c    on the following principles:
c
c          let k = fixed flux level, kp = varying flux level
c          then flux(k)=sum over kp : (deltab(kp)*tau(kp,k))
c               over all kps, from 1 to lp1.
c
c          we can break down the calculations for all ks as follows:
c
c          for all ks k=1 to lp1:
c              flux(k)=sum over kp : (deltab(kp)*tau(kp,k))  (1)
c                      over all kps, from k+1 to lp1
c          and
c              for kp from k+1 to lp1:
c                 flux(kp) = deltab(k)*tau(k,kp)              (2)
c
c          now if tau(k,kp)=tau(kp,k) (symmetrical arrays)
c          we can compute a 1-dimensional array tau1d(kp) from
c          k+1 to lp1, each time k is incremented.
c          equations (1) and (2) then become:
c
c             tau1d(kp) = (values for tau(kp,k) at the particular k)
c             flux(k) = sum over kp : (deltab(kp)*tau1d(kp))   (3)
c             flux(kp) = deltab(k)*tau1d(kp)                   (4)
c
c         the terms for tau (k,k) and other special terms (for
c         nearby layers) must, of course, be handled separately, and
c         with care.
c
c      compute "upper triangle" transmission functions for
c      the 9.6 um band (to3sp) and the 15 um band (over1d). also,
c      the
c      stage 1...compute o3 ,over transmission fctns and avephi
c---do k=1 calculation (from flux layer kk to the top) separately
c   as vectorization is improved,and ozone cts transmissivity
c   may be extracted here.
      do 3021 k=1,l
      do 3021 i=1,ipts
      avephi(i,k)=totphi(i,k+1)
3021  continue
c---in order to properly evaluate emiss integrated over the (lp1)
c   layer, a special evaluation of emiss is done. this requires
c   a special computation of avephi, and it is stored in the
c   (otherwise vacant) lp1th position
c
      do 803 i=1,ipts
      avephi(i,lp1)=avephi(i,lm1)+emx1(i)
803   continue
c   compute fluxes for k=1
      call e1e290(ipts,e1cts1,e1cts2,e1flx,e1ctw1,e1ctw2,emiss,
     1            fxo,dt,fxoe2,dte2,avephi,temp,t)
      do 302 k=1,l
      do 302 i=1,ipts
      fac1(i,k)=bo3rnd(2)*tphio3(i,k+1)/toto3(i,k+1)
      to3spc(i,k)=haf*(fac1(i,k)*
     1    (sqrt(one+(four*ao3rnd(2)*toto3(i,k+1))/fac1(i,k))-one))
c   for k=1, to3sp is used instead of to31d (they are equal in this
c   case); to3sp is passed to spa90, while to31d is a work-array.
      to3sp(i,k)=exp(hm1ez*(to3spc(i,k)+sko3r*totvo2(i,k+1)))
      over1d(i,k)=exp(hm1ez*(sqrt(ab15wd*totphi(i,k+1))+
     1            skc1r*totvo2(i,k+1)))
c---because all continuum transmissivities are obtained from the
c  2-d quantity cnttau (and its reciprocal totevv) we store both
c  of these here. for k=1, cont1d equals cnttau
      cnttau(i,k)=exp(hm1ez*totvo2(i,k+1))
      totevv(i,k)=1./cnttau(i,k)
302   continue
      do 3022 k=1,l
      do 3022 i=1,ipts
      co2sp(i,k+1)=over1d(i,k)*co21(i,1,k+1)
3022  continue
      do 3023 k=1,l
      do 3023 i=1,ipts
      co21(i,k+1,1)=co21(i,k+1,1)*over1d(i,k)
3023  continue
c---rlog is the nbl amount for the 15 um band calculation
      do 1808 i=1,ipts
      rlog(i,1)=over1d(i,1)*co2nbl(i,1)
1808  continue
c---the terms when kp=1 for all k are the photon exchange with
c   the top of the atmosphere, and are obtained differently than
c   the other calculations
      do 305 k=2,lp1
      do 305 i=1,ipts
      flx(i,k)= (tc(i,1)*e1flx(i,k)
     1          +ss1(i,1)*cnttau(i,k-1)
     2          +sorc(i,1,13)*to3sp(i,k-1)
     3          +csour(i,1)*co2sp(i,k))
     4          *cldfac(i,1,k)
305   continue
      do 307 i=1,ipts
      flx(i,1)= tc(i,1)*e1flx(i,1)+ss1(i,1)+sorc(i,1,13)
     1          +csour(i,1)
307   continue
c---the kp terms for k=1...
      do 303 kp=2,lp1
      do 303 i=1,ipts
      flx(i,1)=flx(i,1)+(oss(i,kp)*to3sp(i,kp-1)
     1                  +ss2(i,kp)*cnttau(i,kp-1)
     2                  +css(i,kp)*co21(i,kp,1)
     3                  +dtc(i,kp)*emiss(i,kp-1))*cldfac(i,kp,1)
303   continue
c...    ditto for clear sky....
#ifdef CLR
      do 2305 k=2,lp1
      do 2305 i=1,ipts
      flx0(i,k)=  tc(i,1)*e1flx(i,k)
     1          +ss1(i,1)*cnttau(i,k-1)
     2          +sorc(i,1,13)*to3sp(i,k-1)
     3          +csour(i,1)*co2sp(i,k)
#endif
2305  continue
#ifdef CLR
      do 2307 i=1,ipts
      flx0(i,1)=tc(i,1)*e1flx(i,1)+ss1(i,1)+sorc(i,1,13)
     1          +csour(i,1)
#endif
2307  continue
c---the kp terms for k=1...
#ifdef CLR
      do 2303 kp=2,lp1
       do 2303 i=1,ipts
        flx0(i,1)=flx0(i,1)+ oss(i,kp)*to3sp(i,kp-1)
     1                  +ss2(i,kp)*cnttau(i,kp-1)
     2                  +css(i,kp)*co21(i,kp,1)
     3                  +dtc(i,kp)*emiss(i,kp-1)
#endif
2303  continue
c          subroutine spa88 is called to obtain exact cts for water
c     co2 and o3, and approximate cts co2 and o3 calculations.
c
      call spa88(ipts,excts,ctso3,gxcts,sorc,csour,
#ifdef CLR
     1           excts0,ctso30,gxcts0,
#endif
     1           cldfac,temp,press,var1,var2,
     2           p,delp,delp2,totvo2,to3sp,to3spc,
     3           co2sp1,co2sp2,co2sp)
c
c    this section computes the emissivity cts heating rates for 2
c    emissivity bands: the 0-160,1200-2200 cm-1 band and the 800-
c    990,1070-1200 cm-1 band. the remaining cts comtributions are
c    contained in ctso3, computed in spa88.
c
      do 998 i=1,ipts
      vtmp3(i,1)=1.
998   continue
      do 999 k=1,l
      do 999 i=1,ipts
      vtmp3(i,k+1)=cnttau(i,k)*cldfac(i,k+1,1)
999   continue
      do 1001 k=1,l
      do 1001 i=1,ipts
      cts(i,k)=radcon*delp(i,k)*(tc(i,k)*
     1     (e1ctw2(i,k)*cldfac(i,k+1,1)-e1ctw1(i,k)*cldfac(i,k,1)) +
     2      ss1(i,k)*(vtmp3(i,k+1)-vtmp3(i,k)))
1001  continue
c
      do 1011 k=1,l
      do 1011 i=1,ipts
      vtmp3(i,k)=tc(i,k)*(cldfac(i,k,1)*(e1cts1(i,k)-e1ctw1(i,k)) -
     1                    cldfac(i,k+1,1)*(e1cts2(i,k)-e1ctw2(i,k)))
1011  continue
      do 1012 i=1,ipts
      flx1e1(i)=tc(i,lp1)*cldfac(i,lp1,1)*
     1          (e1cts1(i,lp1)-e1ctw1(i,lp1))
1012  continue
      do 1014 k=1,l
      do 1013 i=1,ipts
      flx1e1(i)=flx1e1(i)+vtmp3(i,k)
1013  continue
1014  continue
c  ... ditto for clear sky  ...
#ifdef CLR
      do 2998 i=1,ipts
       vtmp3(i,1)=1.
#endif
2998  continue
#ifdef CLR
      do 2999 k=1,l
      do 2999 i=1,ipts
       vtmp3(i,k+1)=cnttau(i,k)
#endif
2999  continue
#ifdef CLR
      do 2001 k=1,l
      do 2001 i=1,ipts
       cts0(i,k)=radcon*delp(i,k)*(tc(i,k)*
     1     (e1ctw2(i,k)-e1ctw1(i,k)) +
     2      ss1(i,k)*(vtmp3(i,k+1)-vtmp3(i,k)))
#endif
2001  continue
c
#ifdef CLR
      do 2011 k=1,l
      do 2011 i=1,ipts
       vtmp3(i,k)=tc(i,k)*(e1cts1(i,k)-e1ctw1(i,k)  -
     1                    (e1cts2(i,k)-e1ctw2(i,k)))
#endif
2011  continue
#ifdef CLR
      do 2012 i=1,ipts
       flx1e0(i)=tc(i,lp1)*
     1          (e1cts1(i,lp1)-e1ctw1(i,lp1))
#endif
2012  continue
#ifdef CLR
      do 2014 k=1,l
       do 2013 i=1,ipts
        flx1e0(i)=flx1e0(i)+vtmp3(i,k)
#endif
2013   continue
2014  continue
c
c---now repeat flux calculations for the k=2..lm1  cases.
c   calculations for flux level l and lp1 are done separately, as all
c   emissivity and co2 calculations are special cases or nearby layers.
c
      do 321 k=2,lm1
      klen=k
c
      do 3218 kk=1,lp1-k
      do 3218 i=1,ipts
      avephi(i,kk+k-1)=totphi(i,kk+k)-totphi(i,k)
3218  continue
      do 1803 i=1,ipts
      avephi(i,lp1)=avephi(i,lm1)+emx1(i)
1803   continue
c---compute emissivity fluxes (e2) for this case. note that
c   we have omitted the nearby later case (emiss(i,k,k)) as well
c   as all cases with k=l or lp1. but these cases have always
c   been handled as special cases, so we may as well compute
c    their fluxes separastely.
c
      call e290(ipts,emissb,emiss,avephi,klen,fxoe2,dte2)
      do 322 kk=1,lp1-k
      do 322 i=1,ipts
      avmo3(i,kk+k-1)=toto3(i,kk+k)-toto3(i,k)
      avpho3(i,kk+k-1)=tphio3(i,kk+k)-tphio3(i,k)
      avvo2(i,kk+k-1)=totvo2(i,kk+k)-totvo2(i,k)
      cont1d(i,kk+k-1)=cnttau(i,kk+k-1)*totevv(i,k-1)
322   continue
c
      do 3221 kk=1,lp1-k
      do 3221 i=1,ipts
      fac1(i,kk+k-1)=bo3rnd(2)*avpho3(i,kk+k-1)/avmo3(i,kk+k-1)
      vtmp3(i,kk+k-1)=haf*(fac1(i,kk+k-1)*
     1  (sqrt(one+(four*ao3rnd(2)*avmo3(i,kk+k-1))/fac1(i,kk+k-1))-one))
      to31d(i,kk+k-1)=exp(hm1ez*(vtmp3(i,kk+k-1)+sko3r*avvo2(i,kk+k-1)))
      over1d(i,kk+k-1)=exp(hm1ez*(sqrt(ab15wd*avephi(i,kk+k-1))+
     1            skc1r*avvo2(i,kk+k-1)))
      co21(i,kk+k,k)=over1d(i,kk+k-1)*co21(i,kk+k,k)
3221  continue
      do 3223 kp=k+1,lp1
      do 3223 i=1,ipts
      co21(i,k,kp)=over1d(i,kp-1)*co21(i,k,kp)
3223  continue
c---rlog is the nbl amount for the 15 um band calculation
      do 1804 i=1,ipts
      rlog(i,k)=over1d(i,k)*co2nbl(i,k)
1804  continue
c---the kp terms for arbirrary k..
      do 3423 kp=k+1,lp1
      do 3423 i=1,ipts
      flx(i,k)=flx(i,k)+(oss(i,kp)*to31d(i,kp-1)
     1                  +ss2(i,kp)*cont1d(i,kp-1)
     2                  +css(i,kp)*co21(i,kp,k)
     3                  +dtc(i,kp)*emiss(i,kp-1))*cldfac(i,kp,k)
3423  continue
      do 3425 kp=k+1,lp1
      do 3425 i=1,ipts
      flx(i,kp)=flx(i,kp)+(oss(i,k)*to31d(i,kp-1)
     1                   +ss2(i,k)*cont1d(i,kp-1)
     2                   +css(i,k)*co21(i,k,kp)
     3                   +dtc(i,k)*emissb(i,kp-1))*cldfac(i,k,kp)
3425  continue
c ....   ditto for clear sky .. cldfac=1.
#ifdef CLR
      do 2423 kp=k+1,lp1
       do 2423 i=1,ipts
        flx0(i,k)=flx0(i,k)+ oss(i,kp)*to31d(i,kp-1)
     1                  +ss2(i,kp)*cont1d(i,kp-1)
     2                  +css(i,kp)*co21(i,kp,k)
     3                  +dtc(i,kp)*emiss(i,kp-1)
#endif
2423  continue
#ifdef CLR
      do 2425 kp=k+1,lp1
       do 2425 i=1,ipts
        flx0(i,kp)=flx0(i,kp)+ oss(i,k)*to31d(i,kp-1)
     1                   +ss2(i,k)*cont1d(i,kp-1)
     2                   +css(i,k)*co21(i,k,kp)
     3                   +dtc(i,k)*emissb(i,kp-1)
#endif
2425  continue
321   continue
c
c   now do k=l case. since the kp loop is length 1, many simplifi-
c   cations occur. also, the co2 quantities (as well as the emiss
c  quantities) are computed in the nbl sedction; therefore, we want
c  only over,to3 and cont1d (over(i,l),to31d(i,l) and cont1d(i,l)
c  according to the notation. thus no call is made to the e290
c  subroutine.
c         the third section calculates boundary layer and nearby layer
c     corrections to the transmission functions obtained above. methods
c     are given in ref. (4).
c          the following ratios are used in various nbl calculations:
c
c   the remaining calculations are for :
c                        1) the (k,k) terms, k=2,lm1;
c                        2) the (l,l) term
c                        3) the (l,lp1) term
c                        4) the (lp1,l) term
c                        5) the (lp1,lp1) term.
c     each is uniquely handled; different flux terms are computed
c     differently
c
c
c          fourth section obtains water transmission functions
c     used in q(approx) calculations and also makes nbl corrections:
c     1) emiss (i,j) is the transmission function matrix obtained
c     by calling subroutine e1e288;
c     2) "nearby layer" corrections (emiss(i,i)) are obtained
c     using subroutine e3v88;
c     3) special values at the surface (emiss(l,lp1),emiss(lp1,l),
c     emiss(lp1,lp1)) are calculated.
c
c
c      obtain arguments for e1e288 and e3v88:
c
      do 821 i=1,ipts
      tpl(i,1)=temp(i,l)
      tpl(i,lp1)=haf*(t(i,lp1)+temp(i,l))
      tpl(i,llp1)=haf*(t(i,l)+temp(i,l))
821   continue
      do 823 k=2,l
      do 823 i=1,ipts
      tpl(i,k)=t(i,k)
      tpl(i,k+l)=t(i,k)
823   continue
c
c---e2 functions are required in the nbl calculations for 2 cases,
c   denoted (in old code) as (l,lp1) and (lp1,lp1)
      do 833 i=1,ipts
      avephi(i,1)=var2(i,l)
      avephi(i,2)=var2(i,l)+empl(i,l)
833   continue
      call e2spec(ipts,emiss,avephi,fxosp,dtsp)
c
c     call e3v88 for nbl h2o transmissivities
           call e3v88(ipts,emd,tpl,empl)
c
c   compute nearby layer and special-case transmissivities for emiss
c    using methods for h2o given in ref. (4)
      do 851 k=2,l
      do 851 i=1,ipts
      emisdg(i,k)=emd(i,k+l)+emd(i,k)
851   continue
c
c   note that emx1/2 (pressure scaled paths) are now computed in
c   lwr88
      do 861 i=1,ipts
      emspec(i,1)=(emd(i,1)*empl(i,1)-emd(i,lp1)*empl(i,lp1))/
     1 emx1(i) + quartr*(emiss(i,1)+emiss(i,2))
      emisdg(i,lp1)=two*emd(i,lp1)
      emspec(i,2)=two*(emd(i,1)*empl(i,1)-emd(i,llp1)*empl(i,llp1))/
     * emx2(i)
861   continue
      do 331 i=1,ipts
      fac1(i,l)=bo3rnd(2)*var4(i,l)/var3(i,l)
      vtmp3(i,l)=haf*(fac1(i,l)*
     1    (sqrt(one+(four*ao3rnd(2)*var3(i,l))/fac1(i,l))-one))
      to31d(i,l)=exp(hm1ez*(vtmp3(i,l)+sko3r*cntval(i,l)))
      over1d(i,l)=exp(hm1ez*(sqrt(ab15wd*var2(i,l))+
     1            skc1r*cntval(i,l)))
      cont1d(i,l)=cnttau(i,l)*totevv(i,lm1)
      rlog(i,l)=over1d(i,l)*co2nbl(i,l)
331   continue
      do 618 k=1,l
      do 618 i=1,ipts
      rlog(i,k)=log(rlog(i,k))
618   continue
      do 601 k=1,lm1
      do 601 i=1,ipts
      delpr1(i,k+1)=delp(i,k+1)*(press(i,k+1)-p(i,k+1))
      alp(i,k+l)=-sqrt(delpr1(i,k+1))*rlog(i,k+1)
601   continue
      do 603 k=1,l
      do 603 i=1,ipts
      delpr2(i,k+1)=delp(i,k)*(p(i,k+1)-press(i,k))
      alp(i,k)=-sqrt(delpr2(i,k+1))*rlog(i,k)
603   continue
      do 625 i=1,ipts
      alp(i,ll)=-rlog(i,l)
      alp(i,llp1)=-rlog(i,l)*sqrt(delp(i,l)*(p(i,lp1)-press(i,lm1)))
625   continue
c        the first computation is for the 15 um band,with the
c     for the combined h2o and co2 transmission function.
c
c       perform nbl computations for the 15 um band
c***the statement function sf in prev. versions is now explicitly
c   evaluated.
      do 631 k=1,llp1
      do 631 i=1,ipts
      c(i,k)=alp(i,k)*(hmp66667+alp(i,k)*(quartr+alp(i,k)*hm6666m2))
631   continue
      do 641 i=1,ipts
      co21(i,lp1,lp1)=one+c(i,l)
      co21(i,lp1,l)=one+(delp2(i,l)*c(i,ll)-(press(i,l)-p(i,l))*
     1 c(i,llm1))/(p(i,lp1)-press(i,l))
      co21(i,l,lp1)=one+((p(i,lp1)-press(i,lm1))*c(i,llp1)-
     1 (p(i,lp1)-press(i,l))*c(i,l))/(press(i,l)-press(i,lm1))
641   continue
      do 643 k=2,l
      do 643 i=1,ipts
      co21(i,k,k)=one+haf*(c(i,lm1+k)+c(i,k-1))
643   continue
c
c    compute nearby-layer transmissivities for the o3 band and for the
c    one-band continuum band (to3 and emiss2). the sf2 function is
c    used. the method is the same as described for co2 in ref (4).
      do 651 k=1,lm1
      do 651 i=1,ipts
      csub(i,k+1)=cntval(i,k+1)*delpr1(i,k+1)
      csub(i,k+l)=cntval(i,k)*delpr2(i,k+1)
651   continue
c---the sf2 function in prev. versions is now explicitly evaluated
      do 655 k=1,llm2
      do 655 i=1,ipts
      csub2(i,k+1)=sko3r*csub(i,k+1)
      c(i,k+1)=csub(i,k+1)*(hmp5+csub(i,k+1)*
     1          (hp166666-csub(i,k+1)*h41666m2))
      c2(i,k+1)=csub2(i,k+1)*(hmp5+csub2(i,k+1)*
     1           (hp166666-csub2(i,k+1)*h41666m2))
655   continue
      do 661 i=1,ipts
      contdg(i,lp1)=1.+c(i,llm1)
      to3dg(i,lp1)=1.+c2(i,llm1)
661   continue
      do 663 k=2,l
      do 663 i=1,ipts
      contdg(i,k)=one+haf*(c(i,k)+c(i,lm1+k))
      to3dg(i,k)=one+haf*(c2(i,k)+c2(i,lm1+k))
663   continue
c---now obtain fluxes
c
c    for the diagonal terms...
      do 871 k=2,lp1
      do 871 i=1,ipts
      flx(i,k)=flx(i,k)+(dtc(i,k)*emisdg(i,k)
     1                 +ss2(i,k)*contdg(i,k)
     2                 +oss(i,k)*to3dg(i,k)
     3                 +css(i,k)*co21(i,k,k))*cldfac(i,k,k)
871   continue
c     for the two off-diagonal terms...
      do 873 i=1,ipts
      flx(i,l)=flx(i,l)+(css(i,lp1)*co21(i,lp1,l)
     1                  +dtc(i,lp1)*emspec(i,2)
     2                  +oss(i,lp1)*to31d(i,l)
     3                  +ss2(i,lp1)*cont1d(i,l))*cldfac(i,lp1,l)
      flx(i,lp1)=flx(i,lp1)+(css(i,l)*co21(i,l,lp1)
     1                      +oss(i,l)*to31d(i,l)
     2                      +ss2(i,l)*cont1d(i,l)
     3                      +dtc(i,l)*emspec(i,1))*cldfac(i,l,lp1)
873   continue
c ...   ditto for clear sky ... cldfac =1.
#ifdef CLR
      do 2871 k=2,lp1
       do 2871 i=1,ipts
        flx0(i,k)=flx0(i,k)+ dtc(i,k)*emisdg(i,k)
     1                 +ss2(i,k)*contdg(i,k)
     2                 +oss(i,k)*to3dg(i,k)
     3                 +css(i,k)*co21(i,k,k)
#endif
2871  continue
c     for the two off-diagonal terms...
#ifdef CLR
      do 2873 i=1,ipts
       flx0(i,l)=flx0(i,l)+ css(i,lp1)*co21(i,lp1,l)
     1                  +dtc(i,lp1)*emspec(i,2)
     2                  +oss(i,lp1)*to31d(i,l)
     3                  +ss2(i,lp1)*cont1d(i,l)
       flx0(i,lp1)=flx0(i,lp1)+ css(i,l)*co21(i,l,lp1)
     1                      +oss(i,l)*to31d(i,l)
     2                      +ss2(i,l)*cont1d(i,l)
     3                      +dtc(i,l)*emspec(i,1)
#endif
2873  continue
c
c     final section obtains emissivity heating rates,
c     total heating rates and the flux at the ground
c
c     .....calculate the emissivity heating rates
      do 1101 k=1,l
      do 1101 i=1,ipts
      heatem(i,k)=radcon*(flx(i,k+1)-flx(i,k))*delp(i,k)
1101  continue
c     .....calculate the total heating rates
      do 1103 k=1,l
      do 1103 i=1,ipts
      heatra(i,k)=heatem(i,k)-cts(i,k)-ctso3(i,k)+excts(i,k)
1103  continue
c     .....calculate the flux at each flux level using the flux at the
c    top (flx1e1+gxcts) and the integral of the heating rates (vsum1)
      do 1111 k=1,l
      do 1111 i=1,ipts
      vsum1(i,k)=heatra(i,k)*delp2(i,k)*radcon1
1111  continue
      do 1115 i=1,ipts
      topflx(i)=flx1e1(i)+gxcts(i)
      flxnet(i,1)=topflx(i)
1115  continue
c---only the surface value of flux (grnflx) is needed unless
c    the thick cloud section is invoked.
      do 1123 k=2,lp1
      do 1123 i=1,ipts
      flxnet(i,k)=flxnet(i,k-1)+vsum1(i,k-1)
1123  continue
      do 1125 i=1,ipts
      grnflx(i)=flxnet(i,lp1)
1125  continue
c ...   ditto for clear sky .. cldfac=1.
#ifdef CLR
      do 2101 k=1,l
      do 2101 i=1,ipts
       heatem(i,k)=radcon*(flx0(i,k+1)-flx0(i,k))*delp(i,k)
#endif
2101  continue
c     .....calculate the total heating rates
#ifdef CLR
      do 2103 k=1,l
      do 2103 i=1,ipts
       heatr0(i,k)=heatem(i,k)-cts0(i,k)-ctso30(i,k)+excts0(i,k)
#endif
2103  continue
c     .....calculate the flux at each flux level using the flux at the
c    top (flx1e1+gxcts) and the integral of the heating rates (vsum1)
#ifdef CLR
      do 2111 k=1,l
      do 2111 i=1,ipts
       vsum1(i,k)=heatr0(i,k)*delp2(i,k)*radcon1
#endif
2111  continue
#ifdef CLR
      do 2115 i=1,ipts
       topfx0(i)=flx1e0(i)+gxcts0(i)
       flxnt0(i,1)=topfx0(i)
#endif
2115  continue
c---only the surface value of flux (grnflx) is needed unless
c    the thick cloud section is invoked.
#ifdef CLR
      do 2123 k=2,lp1
       do 2123 i=1,ipts
        flxnt0(i,k)=flxnt0(i,k-1)+vsum1(i,k-1)
#endif
2123  continue
#ifdef CLR
      do 2125 i=1,ipts
       grnfx0(i)=flxnt0(i,lp1)
#endif
2125  continue
c
c     this is the thick cloud section.optionally,if thick cloud
c     fluxes are to be "convectively adjusted",ie,df/dp is constant,
c     for cloudy part of grid point, the following code is executed.
c***first,count the number of clouds along the lat. row. skip the
c   entire thick cloud computation of there are no clouds.
c     icnt=0
c     do 1301 i=1,ipts
c     icnt=icnt+nclds(i)
c301  continue
c     if (icnt.eq.0) go to 6999
c---find the maximum number of clouds in the latitude row
c     kclds=nclds(1)
c     do 2106 i=2,ipts
c     kclds=max(nclds(i),kclds)
c106  continue
c
c
c***obtain the pressures and fluxes of the top and bottom of
c   the ncth cloud (it is assumed that all ktop and kbtms have
c   been defined!).
c     do 1361 kk=1,kclds
c     kmin=lp1
c     kmax=0
c     do 1362 i=1,ipts
c       j1=ktop(i,kk+1)
c       if (j1.eq.1) go to 1362
c       j3=kbtm(i,kk+1)
c       if (j3.gt.j1) then
c         ptop(i)=p(i,j1)
c         pbot(i)=p(i,j3+1)
c         ftop(i)=flxnet(i,j1)
c         fbot(i)=flxnet(i,j3+1)
c***obtain the "flux derivative" df/dp (delptc)
c         delptc(i)=(ftop(i)-fbot(i))/(ptop(i)-pbot(i))
c         kmin=min(kmin,j1)
c         kmax=max(kmax,j3)
c       endif
c362  continue
c     kmin=kmin+1
c***calculate the tot. flux chg. from the top of the cloud, for
c   all levels.
c     do 1365 k=kmin,kmax
c     do 1363 i=1,ipts
c       if (ktop(i,kk+1).eq.1) go to 1363
c       if(ktop(i,kk+1).lt.k .and. k.le.kbtm(i,kk+1)) then
c         z1(i,k)=(p(i,k)-ptop(i))*delptc(i)+ftop(i)
coriginal flxnet(i,k)=flxnet(i,k)*(one-camt(i,kk+1)) +
coriginal1            z1(i,k)*camt(i,kk+1)
c         flxnet(i,k)=z1(i,k)
c       endif
c363  continue
c365  continue
c361  continue
c***using this flux chg. in the cloudy part of the grid box, obtain
c   the new fluxes, weighting the clear and cloudy fluxes:again, only
c    the fluxes in thick-cloud levels will eventually be used.
c     do 6051 k=1,lp1
c     do 6051 i=1,ipts
c     flxnet(i,k)=flxnet(i,k)*(one-camt(i,nc)) +
c    1            z1(i,k)*camt(i,nc)
c051  continue
c***merge flxthk into flxnet for appropriate levels.
c     do 1401 k=1,lp1
c     do 1401 i=1,ipts
c     if (k.gt.itop(i) .and. k.le.ibot(i)
c    1  .and.  (nc-1).le.nclds(i))  then
c          flxnet(i,k)=flxthk(i,k)
c     endif
c401  continue
c
c******end of cloud loop*****
c6001  continue
c6999  continue
c***the final step is to recompute the heating rates based on the
c   revised fluxes:
c      do 6101 k=1,l
c      do 6101 i=1,ipts
c      heatra(i,k)=radcon*(flxnet(i,k+1)-flxnet(i,k))*delp(i,k)
c6101  continue
c     the thick cloud section ends here.
      return
      end
