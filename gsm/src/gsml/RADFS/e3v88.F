#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine e3v88(ipts,emv,tv,av)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
#include <paramodel.h>
#include <hcon.h>
#include <rdparm.h>
#include <tabcom.h>
      dimension it(imbx,llp1),ww1(imbx,llp1),
     1          dt(imbx,llp1),ww2(imbx,llp1),
     2          du(imbx,llp1)
c   the following arrays are equivalenced to vtemp arrays
      dimension fxo(imbx,llp1),fyo(imbx,llp1),tmp3(imbx,llp1)
c
c    dimensions of arrays in argument list
      dimension emv(imbx,llp1),tv(imbx,llp1),av(imbx,llp1)
c
c   the following array is equivalenced to an array in tabcom.h
      dimension em3v(5040)
c
      equivalence (em3v(1),em3(1,1))
      equivalence (fxo,ww1),(fyo,ww2),(it,tmp3)
c---the following loop replaces a double loop over i (1-imax) and
c   k (1-llp1)
      do 203 k=1,llp1
      do 203 i=1,ipts
        fxo(i,k)=aint(tv(i,k)*hp1)
        tmp3(i,k)=log10(av(i,k))+h16e1
        dt(i,k)=tv(i,k)-ten*fxo(i,k)
        fyo(i,k)=aint(tmp3(i,k)*ten)
        du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
c---obtain index for table lookup; this value will have to be
c   decremented by 9 to account for table temps starting at 100k.
        it(i,k)=fxo(i,k)+fyo(i,k)*h28e1
        ww1(i,k)=ten-dt(i,k)
        ww2(i,k)=hp1-du(i,k)
        emv(i,k)=ww1(i,k)*ww2(i,k)*em3v(it(i,k)-9)+
     1           ww2(i,k)*dt(i,k)*em3v(it(i,k)-8)+
     2           ww1(i,k)*du(i,k)*em3v(it(i,k)+19)+
     3           dt(i,k)*du(i,k)*em3v(it(i,k)+20)
203   continue
      return
      end
