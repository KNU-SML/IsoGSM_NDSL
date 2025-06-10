#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine cnvcldn(ims2,imx2,clstp,
     1                  rn,kbot,ktop,cv,cvb,cvt,count)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c-paramodel
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    cnvcld      computes convective cloud cover
c   prgmmr: fanglin_yang     org: cmb/emc/ncep date 00-05-02
c   prgmmr: yifeng cui mp    org: sdsc         date:02-09-05
c
c
c abstract: computes convective cloud cover and cloud tops and bottoms
c   after the deep convection is invoked. cloud cover is based on
c   the emperical relation between cloud cover and  precipitation rate
c   from eq(4) of slingo and slingo (1991, jgr, 15341-15357).
c   [see also slingo (1987), qjrms, 899-927]
c
c   input argument list:
c     clstp    - if >0  interval in hour of accumulation
c                if <=0 flag for initialization
c     rn       - convective rain in meters
c     kbot     - cloud bottom level
c     ktop     - cloud top level
c
c   output argument list:
c     cv       - convective cloud cover (0~1.0)
c     cvb      - convective cloud base level
c     cvt      - convective cloud top level
c
c$$$
      dimension rn(imx2),kbot(imx2),ktop(imx2),
     1          cv(imx2),cvb(imx2),cvt(imx2)
c  local work variables and arrays
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      data cvb0/100./
c
      a0    = 0.200    !coefficient of eq(4) from slingos(1991) (default 0.245)
      b0    = 0.120    !coefficient of eq(4) from slingos(1991) (default 0.125)
      cvmax = 0.75     !upper limit of convectice cloud cover   (default 0.8)
c
      im=ims2
c-----------------------------------------------------------------------
c  initialize convective rain and range
      if(clstp.le.0..and.clstp.gt.-10.) then
        do i=1,im
          cv(i)=0.
          cvb(i)=cvb0
          cvt(i)=0.
        enddo
      endif
c-----------------------------------------------------------------------
c  accumulate convective rain and range
      if(clstp.gt.-99. .and. count.ne.0. ) then
        do i=1,im
          if(rn(i).gt.0.) then
            cv(i)=cv(i)+rn(i)
            cvb(i)=min(cvb(i),float(kbot(i)))
            cvt(i)=max(cvt(i),float(ktop(i)))
          endif
        enddo
      endif
c-----------------------------------------------------------------------
c  convert precipitation rate into cloud fraction
      if (clstp.gt.0..or.(clstp.lt.0.and.clstp.gt.-10.)) then
        do i=1,im
          cvtmp=cv(i)
          if(cvtmp.le.0.) then
            cvb(i)=cvb0
            cvt(i)=0.
            cv(i) =0.
          else
            xx=cvtmp*24000./abs(clstp)  !precip from meter to mm/day
            cv(i)=a0+b0*alog(xx)
            cv(i)=max(0.,min(cv(i),cvmax))
          endif
        enddo
      endif
c-----------------------------------------------------------------------
      return
      end
