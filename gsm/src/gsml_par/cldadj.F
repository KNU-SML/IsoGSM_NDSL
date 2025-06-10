#define PI pi
#undef ADJCNV
#include <define.h>
      subroutine cldadj(cldary,imx22,imx2,kmx,ptop1,ilots,si,xlatrd)
c
c  adjust full level of model cloudiness using latitude dependent ratio table
c
      dimension cldary(imx22,kmx)
      dimension ptop1(ilots,4)
      dimension si(kmx+1)
      dimension xlatrd(imx2)
c
      parameter(idmr=19)
      dimension ratio(idmr,3)
c
      logical lstr
c
c     data ratio/1.2, 1.4, 1.4, 1.2, 1.1, 1.2, 1.8, 2.1, 1.8, 1.3,
c    1           1.5, 1.7, 1.5, 1.3, 1.1, 1.0, 1.0, 1.1, 1.3,
c    2           0.8, 0.7, 0.8, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6,
c    3           0.5, 0.5, 0.7, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5,
c    4           0.4, 0.2, 0.4, 0.4, 0.5, 0.7, 0.6, 0.5, 0.7, 0.7,
c    5           0.5, 0.4, 0.6, 0.8, 0.6, 0.2, 0.2, 0.2, 0.2/
c
c 1st iteration
c     data ratio/1.0, 0.2, 0.6, 0.6, 0.6, 1.0, 1.5, 1.4, 1.1, 0.6,
c    1           0.8, 1.3, 1.4, 1.0, 0.7, 0.6, 0.8, 0.5, 0.3,
c
c    2           0.1, 0.5, 0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.4, 0.4,
c    3           0.4, 0.3, 0.6, 1.0, 1.1, 1.1, 1.1, 1.3, 1.3,
c
c    4           1.2, 0.9, 0.3, 0.4, 0.6, 0.6, 0.6, 0.5, 0.6, 0.5,
c    5           0.5, 0.3, 0.4, 0.7, 0.6, 0.3, 0.2, 0.2, 0.1/
c final iteration
c     data ratio/1.0, 0.2, 0.6, 0.6, 0.7, 1.5, 1.5, 1.3, 1.1, 0.5,
c    1           1.2, 1.3, 1.4, 1.0, 0.6, 0.5, 0.8, 0.4, 0.2,
c
c    2           0.1, 0.7, 1.3, 0.9, 0.7, 0.6, 0.4, 0.2, 0.2, 0.2,
c    3           0.3, 0.2, 0.5, 1.0, 1.2, 1.1, 1.2, 1.4, 1.1,
c
c    4           2.0, 1.6, 0.1, 0.2, 0.5, 0.6, 0.5, 0.3, 0.5, 0.3,
c    5           0.5, 0.2, 0.6, 0.6, 0.4, 0.1, 0.1, 0.1, 0.1/
c yet final iteration
c     data ratio/1.5, 1.5, 1.5, 0.6, 0.7, 1.8, 1.3, 1.4, 1.0, 0.5,
c    1           1.2, 1.2, 1.3, 1.0, 0.5, 0.6, 0.8, 0.8, 0.8,
c
c    2           1.6, 1.6, 1.6, 0.8, 0.7, 0.5, 0.3, 0.1, 0.1, 0.1,
c    3           0.3, 0.1, 0.4, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2,
c
c    4           0.05, 0.05, 0.05, 0.1, 0.6, 0.7, 0.6, 0.3, 0.6, 0.4,
c    5           0.7, 0.2, 1.0, 0.6, 0.3, 0.1, 0.05, 0.05, 0.05/
c
      data ratio/0.4, 0.4, 0.6, 0.6, 0.7, 1.5, 1.5, 1.4, 1.0, 1.0,
     1           1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.4, 0.2,
 
     2           0.5, 0.5, 0.5, 0.7, 0.7, 0.5, 0.3, 0.02, 0.02, 0.01,
     3           0.01, 0.02, 0.4, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2,
 
     4           0.02, 0.02, 0.02, 0.3, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
     5           0.6, 0.6, 0.6, 0.6, 0.3, 0.05, 0.02, 0.02, 0.02/
c
#include <paramodel.h>
c
      do k=1,kmx 
        do i=1,imx2
c
c  find index of ratio from latitude and level
c
          kc=0
          do kk=1,3
            if(si(k).le.ptop1(i,kk  )*1.0e-3.and.
     1         si(k).gt.ptop1(i,kk+1)*1.e-03) kc=kk
          enddo
          if(kc.gt.0.and.cldary(i,k).gt.0.) then
c
            rlat=float(idmr-1)/pi_*xlatrd(i)+float(idmr+1)/2.
            ilat=int(rlat)
            ilat=min(max(ilat,1),idmr)
c
c interpolate ratio to the model grid
c
            ilatp=min(ilat+1,idmr)
            ratiox=ratio(ilat,kc)+
     1        (ratio(ilatp,kc)-ratio(ilat,kc))*(rlat-float(ilat))
c
c identify and extract convective clouds
c
            if(cldary(i,k).gt.1.e2) then
              ccnv=float(int(cldary(i,k))/10)*1.0e-3
              cstr=cldary(i,k)-10.e0*aint(1.e3*ccnv)
              if(ccnv.lt.0..or.ccnv.gt.1.0) then
                print *,'convective cloud amount incorrect. ccnv=',ccnv
#ifdef MP
#ifdef RSM
                call rmpabort
#else
                call mpabort
#endif
#else
                call abort
#endif
              endif
            else
              ccnv=0.
              cstr=cldary(i,k)
            endif
c
c identify temp inversion clouds
c
            lstr=.false.
            if(cstr.gt.2) then
              lstr=.true.
            endif
            cstr=mod(cstr,2.0e0)
c
            if(cstr.lt.0..or.cstr.gt.1.0) then
              print *,'stable cloud amount incorrect. cstr=',cstr
#ifdef MP
#ifdef RSM
              call rmpabort
#else
              call mpabort
#endif
#else
              call abort
#endif
            endif
c
c  adjust cloud amount
c
#ifdef ADJCNV
#ifdef ICECLOUD
            ccnv=min(max(ccnv*ratiox,0.),1.0)
#else
            ccnv=ccnv*ratiox
#endif
#else
            ccnv=ccnv
#endif
#ifdef ICECLOUD
            cstr=min(max(cstr*ratiox,0.),1.0)
#else
            cstr=cstr*ratiox
#endif
c
            if(lstr) then
              cstr=cstr+2.0e0
            endif
c
c  reconstruct cldary
c
            cldary(i,k)=cstr+10.0e0*aint(1.0e3*ccnv)
c
          endif
        enddo
      enddo
c
      return
      end
