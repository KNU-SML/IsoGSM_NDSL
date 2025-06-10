#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine accdia(lons2,lonf2,a,dt,kd,gda)
c  accumulate diagnostics
c
c  nlonx is 1st dimension of a
c  a is array to accumulate
c  dt is factor by which to multiply before accumulating
c  kd is diagnostic number
c  gda is accumulation array
c
#include <paramodel.h>
#include <comgda.h>
c
#ifdef DG3
#ifndef RSM
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
#else
#ifdef MP
#define LONF2S igrd12p_
#else
#define LONF2S igrd12_
#endif
#endif
      dimension a(lonf2,levs_)
      dimension gda(nwgda,kdgda)
c
      if(kd.gt.0.and.kd.le.kdgda) then
      do 20 k=1,levs_
      do 10 i=1,lons2
      ik=i+LONF2S*(k-1)
      gda(ik,kd)=gda(ik,kd)+a(i,k)*dt
10    continue
20    continue
      endif
#endif
c
      return
      end
