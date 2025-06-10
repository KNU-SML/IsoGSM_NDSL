#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine e1e290(ipts,g1,g2,g3,g4,g5,emiss,fxoe1,dte1,fxoe2,dte2,
     1                  avephi,temp,t)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c     subroutine e1e290 computes the exchange terms in the flux equation
c  for longwave radiation for all terms except the exchange with the
c  top of the atmosphere. the method is a table lookup on a pre-
c  computed e2 function (defined in ref. (4)).
c      the e1 function  calculations (formerly done in subroutine
c  e1v88 compute the flux resulting from the exchange of photons
c  between a layer and the top of the atmosphere.  the method is a
c  table lookup on a pre-computed e1 function.
c     calculations are done in two frequency ranges:
c       1) 0-560,1200-2200 cm-1   for q(approx)
c       2) 160-560 cm-1           for q(approx,cts).
c  motivation for these calculations is in references (1) and (4).
c       inputs:                    (common blocks)
c     table1,table2,table3,em1,em1wde  tabcom
c     avephi                           tfcom
c     temp                             radisw
c     t                                kdacom
c     fxoe1,dte1                argument list
c     fxoe2,dte2                argument list
c       outputs:
c     emiss                            tfcom
c     g1,g2,g3                  argument list,for 1st freq. range
c     g4,g5                     argument list,for 2nd freq. range
c
c        called by :     fst88
c        calls     :
c
c
#include <paramodel.h>
#include <hcon.h>
#include <rdparm.h>
#include <tabcom.h>
      dimension temp(imbx,lp1),t(imbx,lp1)
      dimension avephi(imbx,lp1),emiss(imbx,lp1)
c
      dimension it1(imbx,ll3p),ival(imbx,lp1),
     5  fyo(imbx,lp1),du(imbx,lp1),
     6  ww1(imbx,lp1),ww2(imbx,lp1),
     7  tmp3(imbx,lp1),tmp5(imax),tmp9(imax)
c---variables equivalenced to common block variables
      dimension t1(5040),t2(5040),t4(5040)
      dimension em1v(5040),em1vw(5040)
c---variables in the argument list
      dimension fxoe1(imbx,lp1),dte1(imbx,lp1),
     1          fxoe2(imbx,lp1),dte2(imbx,lp1),
     2    g1(imbx,lp1),g2(imbx,l),g3(imbx,lp1),g4(imbx,lp1),g5(imbx,l)
c
      equivalence (em1v(1),em1(1,1)),(em1vw(1),em1wde(1,1))
      equivalence (t1(1),table1(1,1)),(t2(1),table2(1,1)),
     1 (t4(1),table3(1,1))
c---first we obtain the emissivities as a function of temperature
c   (index fxo) and water amount (index fyo). this part of the code
c   thus generates the e2 function. the fxo indices have been
c   obtained in fst88, for convenience.
c
c---this subroutine evaluates the k=1 case only--
c
c---this loop replaces loops going fromi=1,imax and kp=2,lp1 plus
c   the special case for the lp1th layer.
      do 1322 kp=1,lp1
      do 1322 i=1,ipts
      tmp3(i,kp)=log10(avephi(i,kp))+h16e1
      fyo(i,kp)=aint(tmp3(i,kp)*ten)
      du(i,kp)=tmp3(i,kp)-hp1*fyo(i,kp)
      fyo(i,kp)=h28e1*fyo(i,kp)
      ival(i,kp)=fyo(i,kp)+fxoe2(i,kp)
      emiss(i,kp)=t1(ival(i,kp))+du(i,kp)*t2(ival(i,kp))
     1                        +dte2(i,kp)*t4(ival(i,kp))
1322  continue
c
c---the special case emiss(i,l) (layer kp) is obtained now
c   by averaging the values for l and lp1:
      do 1344 i=1,ipts
      emiss(i,l)=haf*(emiss(i,l)+emiss(i,lp1))
1344  continue
c
c   calculations for the kp=1 layer are not performed, as
c   the radiation code assumes that the top flux layer (above the
c   top data level) is isothermal, and hence contributes nothing
c   to the fluxes at other levels.
c
c***the following is the calculation for the e1 function, formerly
c    done in subroutine e1v88. the move to e1e288 is due to the
c    savings in obtaining index values (the temp. indices have
c    been obtained in fst88, while the u-indices are obtained
c    in the e2 calcs.,with k=1).
c
c
c   for terms involving top layer, du is not known; in fact, we
c   use index 2 to repersent index 1 in prev. code. this means that
c    the it1 index 1 and llp1 has to be calculated separately. the
c   index llp2 gives the same value as 1; it can be omitted.
      do 208 i=1,ipts
      it1(i,1)=fxoe1(i,1)
      ww1(i,1)=ten-dte1(i,1)
      ww2(i,1)=hp1
208   continue
      do 209 kp=1,l
      do 209 i=1,ipts
      it1(i,kp+1)=fyo(i,kp)+fxoe1(i,kp+1)
      it1(i,kp+lp1)=fyo(i,kp)+fxoe1(i,kp)
      ww1(i,kp+1)=ten-dte1(i,kp+1)
      ww2(i,kp+1)=hp1-du(i,kp)
209   continue
      do 211 kp=1,l
      do 211 i=1,ipts
      it1(i,kp+llp1)=fyo(i,kp)+fxoe1(i,1)
211   continue
c
c
c  g3(i,1) has the same values as g1 (and did all along)
      do 230 i=1,ipts
      g1(i,1)=ww1(i,1)*ww2(i,1)*em1v(it1(i,1))+
     1        ww2(i,1)*dte1(i,1)*em1v(it1(i,1)+1)
      g3(i,1)=g1(i,1)
230   continue
      do 240 kp=1,l
      do 240 i=1,ipts
      g1(i,kp+1)=ww1(i,kp+1)*ww2(i,kp+1)*em1v(it1(i,kp+1))+
     1        ww2(i,kp+1)*dte1(i,kp+1)*em1v(it1(i,kp+1)+1)+
     2        ww1(i,kp+1)*du(i,kp)*em1v(it1(i,kp+1)+28)+
     3        dte1(i,kp+1)*du(i,kp)*em1v(it1(i,kp+1)+29)
      g2(i,kp)=ww1(i,kp)*ww2(i,kp+1)*em1v(it1(i,kp+lp1))+
     1        ww2(i,kp+1)*dte1(i,kp)*em1v(it1(i,kp+lp1)+1)+
     1        ww1(i,kp)*du(i,kp)*em1v(it1(i,kp+lp1)+28)+
     2        dte1(i,kp)*du(i,kp)*em1v(it1(i,kp+lp1)+29)
240   continue
      do 241 kp=2,lp1
      do 241 i=1,ipts
      g3(i,kp)=ww1(i,1)*ww2(i,kp)*em1v(it1(i,ll+kp))+
     1        ww2(i,kp)*dte1(i,1)*em1v(it1(i,ll+kp)+1)+
     2        ww1(i,1)*du(i,kp-1)*em1v(it1(i,ll+kp)+28)+
     3        dte1(i,1)*du(i,kp-1)*em1v(it1(i,ll+kp)+29)
241   continue
c
      do 244 i=1,ipts
      g4(i,1)=ww1(i,1)*ww2(i,1)*em1vw(it1(i,1))+
     1        ww2(i,1)*dte1(i,1)*em1vw(it1(i,1)+1)
244   continue
      do 242 kp=1,l
      do 242 i=1,ipts
      g4(i,kp+1)=ww1(i,kp+1)*ww2(i,kp+1)*em1vw(it1(i,kp+1))+
     1        ww2(i,kp+1)*dte1(i,kp+1)*em1vw(it1(i,kp+1)+1)+
     2        ww1(i,kp+1)*du(i,kp)*em1vw(it1(i,kp+1)+28)+
     3        dte1(i,kp+1)*du(i,kp)*em1vw(it1(i,kp+1)+29)
      g5(i,kp)=ww1(i,kp)*ww2(i,kp+1)*em1vw(it1(i,kp+lp1))+
     1        ww2(i,kp+1)*dte1(i,kp)*em1vw(it1(i,kp+lp1)+1)+
     1        ww1(i,kp)*du(i,kp)*em1vw(it1(i,kp+lp1)+28)+
     2        dte1(i,kp)*du(i,kp)*em1vw(it1(i,kp+lp1)+29)
242   continue
c
      return
      end
