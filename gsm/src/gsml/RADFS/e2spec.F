#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine e2spec(ipts,emiss,avephi,fxosp,dtsp)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c     subroutine e2spec computes the exchange terms in the flux equation
c  for longwave radiation for 2 terms used for nearby layer compu-
c  tations. the method is a table lookup on a pre-
c  computed e2 function (defined in ref. (4)).
c     calculations are done in the frequency range:
c        0-560,1200-2200 cm-1
c  motivation for these calculations is in references (1) and (4).
c       inputs:                    (common blocks)
c     table1,table2,table3,            tabcom
c     avephi                           tfcom
c     fxosp,dtsp                argument list
c       outputs:
c     emiss                            tfcom
c
c        called by :     fst88
c        calls     :
c
c
#include <paramodel.h>
#include <hcon.h>
#include <rdparm.h>
#include <tabcom.h>
      dimension avephi(imbx,lp1),emiss(imbx,lp1)
      dimension ival(imbx,lp1),
     1  fyo(imbx,lp1),du(imbx,lp1),
     2  tmp3(imbx,lp1)
c---variables equivalenced to common block variables
      dimension t1(5040),t2(5040),t4(5040)
c---variables in the argument list
      dimension fxosp(imbx,2),dtsp(imbx,2)
c
      equivalence (t1(1),table1(1,1)),(t2(1),table2(1,1)),
     1 (t4(1),table3(1,1))
c---first we obtain the emissivities as a function of temperature
c   (index fxo) and water amount (index fyo). this part of the code
c   thus generates the e2 function.
c
      do 132 k=1,2
      do 132 i=1,ipts
      tmp3(i,k)=log10(avephi(i,k))+h16e1
      fyo(i,k)=aint(tmp3(i,k)*ten)
      du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
      ival(i,k)=h28e1*fyo(i,k)+fxosp(i,k)
      emiss(i,k)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))+
     1                         dtsp(i,k)*t4(ival(i,k))
132   continue
      return
      end
