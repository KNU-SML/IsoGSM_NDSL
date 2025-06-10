#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine e290(ipts,emissb,emiss,avephi,klen,fxoe2,dte2)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c     subroutine e290 computes the exchange terms in the flux equation
c  for longwave radiation for all terms except the exchange with the
c  top of the atmosphere. the method is a table lookup on a pre-
c  computed e2 function (defined in ref. (4)).
c     calculations are done in the frequency range:
c       1) 0-560,1200-2200 cm-1   for q(approx)
c  motivation for these calculations is in references (1) and (4).
c       inputs:                    (common blocks)
c     table1,table2,table3,            tabcom
c     avephi                           tfcom
c     fxoe2,dte2,klen           argument list
c       outputs:
c     emiss,emissb                     tfcom
c
c        called by :     fst88
c        calls     :
c
c
#include <paramodel.h>
#include <hcon.h>
#include <rdparm.h>
#include <tabcom.h>
      dimension emissb(imbx,lp1),emiss(imbx,lp1),avephi(imbx,lp1)
      dimension ival(imbx,lp1),
     1  dt(imbx,lp1),fyo(imbx,lp1),du(imbx,lp1)
c---tmp3 may be equivalenced to dt in vtemp
      dimension tmp3(imbx,lp1)
c---variables equivalenced to common block variables
      dimension t1(5040),t2(5040),t4(5040)
c---variables in the argument list
      dimension fxoe2(imbx,lp1),dte2(imbx,lp1)
c
      equivalence (tmp3,dt)
      equivalence (t1(1),table1(1,1)),(t2(1),table2(1,1)),
     1 (t4(1),table3(1,1))
c---first we obtain the emissivities as a function of temperature
c   (index fxo) and water amount (index fyo). this part of the code
c   thus generates the e2 function.
c
c---calculations for varying kp (from kp=k+1 to lp1, including special
c   case: results are in emiss
      do 132 k=1,lp2-klen
      do 132 i=1,ipts
      tmp3(i,k)=log10(avephi(i,k+klen-1))+h16e1
      fyo(i,k)=aint(tmp3(i,k)*ten)
      du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
      fyo(i,k)=h28e1*fyo(i,k)
      ival(i,k)=fyo(i,k)+fxoe2(i,k+klen-1)
      emiss(i,k+klen-1)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))
     1                           +dte2(i,k+klen-1)*t4(ival(i,k))
132   continue
c---the special case emiss(i,l) (layer kp) is obtained now
c   by averaging the values for l and lp1:
      do 1344 i=1,ipts
      emiss(i,l)=haf*(emiss(i,l)+emiss(i,lp1))
1344  continue
c---note that emiss(i,lp1) is not useful after this point.
c
c---calculations for kp=klen and varying k; results are in emissb.
c  in this case, the temperature index is unchanged, always being
c  fxo(i,klen-1); the water index changes, but is symmetrical with
c  that for the varying kp case.note that the special case is not
c  involved here.
c     (fixed level) k varies from (klen+1) to lp1; results are in
c   emissb(i,(klen) to l)
      do 142 k=1,lp1-klen
      do 142 i=1,ipts
      dt(i,k)=dte2(i,klen-1)
      ival(i,k)=fyo(i,k)+fxoe2(i,klen-1)
142   continue
c
      do 234 k=1,lp1-klen
      do 234 i=1,ipts
      emissb(i,k+klen-1)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))
     1                                +dt(i,k)*t4(ival(i,k))
234   continue
      return
      end
