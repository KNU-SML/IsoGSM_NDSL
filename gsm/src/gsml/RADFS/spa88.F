#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
#include <define.h>
      subroutine spa88(ipts,excts,ctso3,gxcts,sorc,csour,
#ifdef CLR
     1                 excts0,ctso30,gxcts0,
#endif
     1                 cldfac,temp,press,var1,var2,
     2                 p,delp,delp2,totvo2,to3sp,to3spc,
     3                 co2sp1,co2sp2,co2sp)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c-paramodel
c
#include <hcon.h>
#include <rdparm.h>
#include <rnddta.h>
c
      dimension sorc(imbx,lp1,nbly),csour(imbx,lp1)
      dimension cldfac(imbx,lp1,lp1)
      dimension temp(imbx,lp1),press(imbx,lp1)
      dimension var1(imbx,l),var2(imbx,l)
      dimension p(imbx,lp1),delp(imbx,l),delp2(imbx,l)
      dimension totvo2(imbx,lp1),to3spc(imbx,l),to3sp(imbx,lp1)
      dimension co2sp1(imbx,lp1),co2sp2(imbx,lp1),co2sp(imbx,lp1)
      dimension excts(imbx,l),ctso3(imbx,l),gxcts(imax)
#ifdef CLR
      dimension excts0(imbx,l),ctso30(imbx,l),gxcts0(imax)
      dimension ctmp0(imbx,lp1),ctmp20(imbx,lp1),ctmp30(imbx,lp1)
#endif
c
      dimension phitmp(imbx,l),psitmp(imbx,l),tt(imbx,l),
     1          fac1(imbx,l),fac2(imbx,l),
     2          ctmp(imbx,lp1),x(imbx,l),y(imbx,l),
     3          topm(imbx,l),topphi(imbx,l),
     4          ctmp3(imbx,lp1),ctmp2(imbx,lp1)
      dimension f(imbx,l),ff(imbx,l),ag(imbx,l),agg(imbx,l)
      equivalence (f,ag,phitmp)
      equivalence (ff,agg,psitmp)
c---compute temperature quantities for use in program
      do 101 k=1,l
      do 101 i=1,imax
      x(i,k)=temp(i,k)-h25e2
      y(i,k)=x(i,k)*x(i,k)
101   continue
c---initialize ctmp(i,1),ctmp2(i,1),ctmp3(i,1) to unity; these are
c   transmission fctns at the top.
      do 345 i=1,ipts
      ctmp(i,1)=one
      ctmp2(i,1)=1.
      ctmp3(i,1)=1.
345   continue
c...   ditto for the clear sky calculation...
#ifdef CLR
      do 2345 i=1,ipts
       ctmp0(i,1)=one
       ctmp20(i,1)=1.
       ctmp30(i,1)=1.
#endif
2345  continue
c***begin loop on frequency bands (1)***
c
c---calculation for band 1 (combined band 1)
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 301 k=1,l
      do 301 i=1,ipts
      f(i,k)=h44194m2*(apcm(1)*x(i,k)+bpcm(1)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(1)*x(i,k)+btpcm(1)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
301   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 315 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
315   continue
      do 319 k=2,l
      do 317 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
317   continue
319   continue
c---tt is the cloud-free cts transmission function
      do 321 k=1,l
      do 321 i=1,ipts
      fac1(i,k)=acomb(1)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(1)*topphi(i,k))
      tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
321   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 353 k=1,l
      do 353 i=1,ipts
      excts(i,k)=sorc(i,k,1)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=sorc(i,k,1)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
353   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 361 i=1,ipts
      gxcts(i)=                 tt(i,l)*sorc(i,l,1)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,1)-sorc(i,l,1))
#ifdef CLR
       gxcts0(i)=gxcts(i)
#endif
       gxcts(i)=cldfac(i,lp1,1)*gxcts(i)
361   continue
c
c
c-----calculation for band 2 (combined band 2)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 401 k=1,l
      do 401 i=1,ipts
      f(i,k)=h44194m2*(apcm(2)*x(i,k)+bpcm(2)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(2)*x(i,k)+btpcm(2)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
401   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 415 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
415   continue
      do 419 k=2,l
      do 417 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
417   continue
419   continue
c---tt is the cloud-free cts transmission function
      do 421 k=1,l
      do 421 i=1,ipts
      fac1(i,k)=acomb(2)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(2)*topphi(i,k))
      tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
421   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 453 k=1,l
      do 453 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,2)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,2)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
453   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 461 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,2)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,2)-sorc(i,l,2)))
461   continue
#ifdef CLR
      do 2461 i=1,ipts
       gxcts0(i)=gxcts0(i)+               tt(i,l)*sorc(i,l,2)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,2)-sorc(i,l,2))
#endif
2461  continue
c
c-----calculation for band 3 (combined band 3)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 501 k=1,l
      do 501 i=1,ipts
      f(i,k)=h44194m2*(apcm(3)*x(i,k)+bpcm(3)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(3)*x(i,k)+btpcm(3)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
501   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 515 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
515   continue
      do 519 k=2,l
      do 517 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
517   continue
519   continue
c---tt is the cloud-free cts transmission function
      do 521 k=1,l
      do 521 i=1,ipts
      fac1(i,k)=acomb(3)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(3)*topphi(i,k))
      tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
521   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 553 k=1,l
      do 553 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,3)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,3)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
553   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 561 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,3)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,3)-sorc(i,l,3)))
561   continue
#ifdef CLR
      do 2561 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,3)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,3)-sorc(i,l,3))
#endif
2561  continue
c
c-----calculation for band 4 (combined band 4)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 601 k=1,l
      do 601 i=1,ipts
      f(i,k)=h44194m2*(apcm(4)*x(i,k)+bpcm(4)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(4)*x(i,k)+btpcm(4)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
601   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 615 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
615   continue
      do 619 k=2,l
      do 617 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
617   continue
619   continue
c---tt is the cloud-free cts transmission function
      do 621 k=1,l
      do 621 i=1,ipts
      fac1(i,k)=acomb(4)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(4)*topphi(i,k))
      tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
621   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 653 k=1,l
      do 653 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,4)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,4)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
653   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 661 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,4)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,4)-sorc(i,l,4)))
661   continue
#ifdef CLR
      do 2661 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,4)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,4)-sorc(i,l,4))
#endif
2661  continue
c
c-----calculation for band 5 (combined band 5)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 701 k=1,l
      do 701 i=1,ipts
      f(i,k)=h44194m2*(apcm(5)*x(i,k)+bpcm(5)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(5)*x(i,k)+btpcm(5)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
701   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 715 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
715   continue
      do 719 k=2,l
      do 717 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
717   continue
719   continue
c---tt is the cloud-free cts transmission function
      do 721 k=1,l
      do 721 i=1,ipts
      fac1(i,k)=acomb(5)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(5)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(5)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
721   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 753 k=1,l
      do 753 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,5)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,5)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
753   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 761 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,5)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,5)-sorc(i,l,5)))
761   continue
#ifdef CLR
      do 2761 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,5)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,5)-sorc(i,l,5))
#endif
2761  continue
c
c-----calculation for band 6 (combined band 6)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 801 k=1,l
      do 801 i=1,ipts
      f(i,k)=h44194m2*(apcm(6)*x(i,k)+bpcm(6)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(6)*x(i,k)+btpcm(6)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
801   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 815 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
815   continue
      do 819 k=2,l
      do 817 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
817   continue
819   continue
c---tt is the cloud-free cts transmission function
      do 821 k=1,l
      do 821 i=1,ipts
      fac1(i,k)=acomb(6)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(6)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(6)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
821   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 853 k=1,l
      do 853 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,6)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,6)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
853   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 861 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,6)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,6)-sorc(i,l,6)))
861   continue
#ifdef CLR
      do 2861 i=1,ipts
       gxcts0(i)=gxcts0(i)+ tt(i,l)*sorc(i,l,6)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,6)-sorc(i,l,6))
#endif
2861  continue
c
c-----calculation for band 7 (combined band 7)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 901 k=1,l
      do 901 i=1,ipts
      f(i,k)=h44194m2*(apcm(7)*x(i,k)+bpcm(7)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(7)*x(i,k)+btpcm(7)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
901   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 915 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
915   continue
      do 919 k=2,l
      do 917 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
917   continue
919   continue
c---tt is the cloud-free cts transmission function
      do 921 k=1,l
      do 921 i=1,ipts
      fac1(i,k)=acomb(7)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(7)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(7)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
921   continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 953 k=1,l
      do 953 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,7)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,7)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
953   continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 961 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,7)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,7)-sorc(i,l,7)))
961   continue
#ifdef CLR
      do 1961 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,7)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,7)-sorc(i,l,7))
#endif
1961  continue
c
c-----calculation for band 8 (combined band 8)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1001 k=1,l
      do 1001 i=1,ipts
      f(i,k)=h44194m2*(apcm(8)*x(i,k)+bpcm(8)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(8)*x(i,k)+btpcm(8)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1001  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1015 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1015  continue
      do 1019 k=2,l
      do 1017 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1017  continue
1019  continue
c---tt is the cloud-free cts transmission function
      do 1021 k=1,l
      do 1021 i=1,ipts
      fac1(i,k)=acomb(8)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(8)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(8)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1021  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1053 k=1,l
      do 1053 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,8)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,8)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1053  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1061 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,8)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,8)-sorc(i,l,8)))
1061  continue
#ifdef CLR
      do 2061 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,8)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,8)-sorc(i,l,8))
#endif
2061  continue
c
c-----calculation for band 9 ( 560-670 cm-1; includes co2)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1101 k=1,l
      do 1101 i=1,ipts
      f(i,k)=h44194m2*(apcm(9)*x(i,k)+bpcm(9)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(9)*x(i,k)+btpcm(9)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1101  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1115 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1115  continue
      do 1119 k=2,l
      do 1117 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1117  continue
1119  continue
c---tt is the cloud-free cts transmission function
      do 1121 k=1,l
      do 1121 i=1,ipts
      fac1(i,k)=acomb(9)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(9)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(9)*totvo2(i,k+1)*sko2d))*co2sp1(i,k+1)
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1121  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1153 k=1,l
      do 1153 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,9)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,9)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1153  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1161 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,9)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,9)-sorc(i,l,9)))
1161  continue
#ifdef CLR
      do 2161 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,9)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,9)-sorc(i,l,9))
#endif
2161  continue
c
c-----calculation for band 10 (670-800 cm-1; includes co2)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1201 k=1,l
      do 1201 i=1,ipts
      f(i,k)=h44194m2*(apcm(10)*x(i,k)+bpcm(10)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(10)*x(i,k)+btpcm(10)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1201  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1215 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1215  continue
      do 1219 k=2,l
      do 1217 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1217  continue
1219  continue
c---tt is the cloud-free cts transmission function
      do 1221 k=1,l
      do 1221 i=1,ipts
      fac1(i,k)=acomb(10)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(10)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(10)*totvo2(i,k+1)*sko2d))*co2sp2(i,k+1)
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1221  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1253 k=1,l
      do 1253 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,10)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,10)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1253  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1261 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,10)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,10)-sorc(i,l,10)))
1261  continue
#ifdef CLR
      do 3261 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,10)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,10)-sorc(i,l,10))
#endif
3261  continue
c
c-----calculation for band 11 (800-900 cm-1)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1301 k=1,l
      do 1301 i=1,ipts
      f(i,k)=h44194m2*(apcm(11)*x(i,k)+bpcm(11)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(11)*x(i,k)+btpcm(11)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1301  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1315 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1315  continue
      do 1319 k=2,l
      do 1317 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1317  continue
1319  continue
c---tt is the cloud-free cts transmission function
      do 1321 k=1,l
      do 1321 i=1,ipts
      fac1(i,k)=acomb(11)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(11)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(11)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1321  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1353 k=1,l
      do 1353 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,11)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,11)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1353  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1361 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,11)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,11)-sorc(i,l,11)))
1361  continue
#ifdef CLR
      do 3361 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,11)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,11)-sorc(i,l,11))
#endif
3361  continue
c
c-----calculation for band 12 (900-990 cm-1)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1401 k=1,l
      do 1401 i=1,ipts
      f(i,k)=h44194m2*(apcm(12)*x(i,k)+bpcm(12)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(12)*x(i,k)+btpcm(12)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1401  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1415 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1415  continue
      do 1419 k=2,l
      do 1417 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1417  continue
1419  continue
c---tt is the cloud-free cts transmission function
      do 1421 k=1,l
      do 1421 i=1,ipts
      fac1(i,k)=acomb(12)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(12)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(12)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1421  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1453 k=1,l
      do 1453 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,12)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,12)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1453  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1461 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,12)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,12)-sorc(i,l,12)))
1461  continue
#ifdef CLR
      do 3461 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,12)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,12)-sorc(i,l,12))
#endif
3461  continue
c
c-----calculation for band 13 (990-1070 cm-1; includes o3))
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1501 k=1,l
      do 1501 i=1,ipts
      f(i,k)=h44194m2*(apcm(13)*x(i,k)+bpcm(13)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(13)*x(i,k)+btpcm(13)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1501  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1515 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1515  continue
      do 1519 k=2,l
      do 1517 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1517  continue
1519  continue
c---tt is the cloud-free cts transmission function
      do 1521 k=1,l
      do 1521 i=1,ipts
      fac1(i,k)=acomb(13)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(13)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(13)*totvo2(i,k+1)*sko2d +to3spc(i,k)))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1521  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1553 k=1,l
      do 1553 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,13)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,13)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1553  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1561 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,13)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,13)-sorc(i,l,13)))
1561  continue
#ifdef CLR
      do 3561 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,13)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,13)-sorc(i,l,13))
#endif
3561  continue
c
c-----calculation for band 14 (1070-1200 cm-1)
c
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1601 k=1,l
      do 1601 i=1,ipts
      f(i,k)=h44194m2*(apcm(14)*x(i,k)+bpcm(14)*y(i,k))
      ff(i,k)=h44194m2*(atpcm(14)*x(i,k)+btpcm(14)*y(i,k))
      ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
      agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
      phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
      psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
1601  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1615 i=1,ipts
      topm(i,1)=phitmp(i,1)
      topphi(i,1)=psitmp(i,1)
1615  continue
      do 1619 k=2,l
      do 1617 i=1,ipts
      topm(i,k)=topm(i,k-1)+phitmp(i,k)
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
1617  continue
1619  continue
c---tt is the cloud-free cts transmission function
      do 1621 k=1,l
      do 1621 i=1,ipts
      fac1(i,k)=acomb(14)*topm(i,k)
      fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(14)*topphi(i,k))
      tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+
     1           betacm(14)*totvo2(i,k+1)*sko2d))
      ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
1621  continue
c---excts is the cts cooling rate accumulated over frequency bands
      do 1653 k=1,l
      do 1653 i=1,ipts
      excts(i,k)=excts(i,k)+sorc(i,k,14)*
     1             (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,14)*
     1             (ctmp0(i,k+1)-ctmp0(i,k))
#endif
1653  continue
c---gxcts is the exact cts top flux accumulated over frequency bands
      do 1661 i=1,ipts
      gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,14)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,14)-sorc(i,l,14)))
1661  continue
#ifdef CLR
      do 3661 i=1,ipts
       gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,14)+
     1   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +
     2   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     3   (sorc(i,lp1,14)-sorc(i,l,14))
#endif
3661  continue
c
c
c   obtain cts flux at the top by integration of heating rates and
c   using cts flux at the bottom (current value of gxcts). note
c   that the pressure quantities and conversion factors have not
c   been included either in excts or in gxcts. these cancel out, thus
c   reducing computations!
      do 1731 k=1,l
      do 1731 i=1,ipts
      gxcts(i)=gxcts(i)-excts(i,k)
#ifdef CLR
      gxcts0(i)=gxcts0(i)-excts0(i,k)
#endif
1731  continue
c
c   now scale the cooling rate (excts) by including the pressure
c   factor (delp) and the conversion factor (radcon)
      do 1741 k=1,l
      do 1741 i=1,ipts
      excts(i,k)=excts(i,k)*radcon*delp(i,k)
#ifdef CLR
      excts0(i,k)=excts0(i,k)*radcon*delp(i,k)
#endif
1741  continue
c---this is the end of the exact cts computations; at this point
c   excts has its appropriate value.
c
c*** compute approximate cts heating rates for 15um and 9.6 um bands
c     (ctso3)
      do 1711 k=1,l
      do 1711 i=1,ipts
      ctmp2(i,k+1)=co2sp(i,k+1)*cldfac(i,k+1,1)
      ctmp3(i,k+1)=to3sp(i,k)*cldfac(i,k+1,1)
#ifdef CLR
      ctmp20(i,k+1)=co2sp(i,k+1)
      ctmp30(i,k+1)=to3sp(i,k)
#endif
1711  continue
      do 1701 k=1,l
      do 1701 i=1,ipts
      ctso3(i,k)=radcon*delp(i,k)*
     1     (csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) +
     2      sorc(i,k,13)*(ctmp3(i,k+1)-ctmp3(i,k)))
1701  continue
#ifdef CLR
      do 3701 k=1,l
      do 3701 i=1,ipts
       ctso30(i,k)=radcon*delp(i,k)*
     1     (csour(i,k)*(ctmp20(i,k+1)-ctmp20(i,k)) +
     2      sorc(i,k,13)*(ctmp30(i,k+1)-ctmp30(i,k)))
#endif
3701  continue
      return
      end
