#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine cnvcld(ims2,imx2,clstp,
     1                  rn,kbot,ktop,cv,cvb,cvt,count)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c-paramodel
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    cnvcld      computes convective cloud cover
c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
c   prgmmr:  yifeng cui mpi  org:sdsc        date: 09-09-05
c
c
c abstract: computes convective cloud cover and cloud tops and bottoms
c   after the deep convection is invoked. cloud cover is interpolated
c   from a table relating cloud cover to precipitation rate.
c
c program history log:
c   91-05-07  iredell
c
c usage:    call cnvcld(cmean,lat,iistp,dt,rn,kbot,ktop,cv,cvb,cvt)
c
c   input argument list:
c     cmean    - real flag (ge 0 to accumulate, eq 99 to return values)
c     lat      - integer latitude index
c     iistp    - integer time step number
c     dt       - real time step in seconds
c     rn       - real (nx) convective rain in meters
c     kbot     - integer (nx) cloud bottom level
c     ktop     - integer (nx) cloud top level
c
c   output argument list:
c     cv       - real (nx,ny) convective cloud cover
c     cvb      - real (nx,ny) convective cloud base level
c     cvt      - real (nx,ny) convective cloud top level
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
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
      dimension nmd(ILOTS),pmd(ILOTS)
c  local save variables and arrays
      parameter(ncc=9)
      dimension cc(ncc),p(ncc)
      data cc/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/
      data p/.14,.31,.70,1.6,3.4,7.7,17.,38.,85./
      data cvb0/100./
c
      im=ims2
c-----------------------------------------------------------------------
c  initialize convective rain and range
      if(clstp.le.0..and.clstp.gt.-10.) then
        do i=1,im
          cv(i)=0.
          cvb(i)=cvb0
          cvt(i)=0.
c         cvb(i)=0.
c         cvt(i)=0.
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
c           cvt(i)=max(cvt(i),float(ktop(i)+1))
c           cvb(i)=cvb(i)+kbot(i)*rn(i)
c           cvt(i)=cvt(i)+(ktop(i))*rn(i)
c           cvt(i)=cvt(i)+(ktop(i)+1)*rn(i)
          endif
        enddo
      endif
c-----------------------------------------------------------------------
c  convert precipitation rate into cloud fraction
      if (clstp.gt.0..or.(clstp.lt.0.and.clstp.gt.-10.)) then
        do i=1,im
          if(cv(i).gt.0.) then
c           cvb(i)=nint(cvb(i)/cv(i))
c           cvt(i)=nint(cvt(i)/cv(i))
          else
            cvb(i)=cvb0
            cvt(i)=0.
          endif
          pmd(i)=cv(i)*(24.e+3/abs(clstp))
          nmd(i)=0
        enddo
        do n=1,ncc
          do i=1,im
            if(pmd(i).gt.p(n)) nmd(i)=n
          enddo
        enddo
        do i=1,im
          if(nmd(i).eq.0) then
            cv(i)=0.
            cvb(i)=cvb0
            cvt(i)=0.
          elseif(nmd(i).eq.ncc) then
            cv(i)=cc(ncc)
          else
            cc1=cc(nmd(i))
            cc2=cc(nmd(i)+1)
            p1=p(nmd(i))
            p2=p(nmd(i)+1)
            cv(i)=cc1+(cc2-cc1)*(pmd(i)-p1)/(p2-p1)
          endif
        enddo
      endif
c-----------------------------------------------------------------------
      return
      end
