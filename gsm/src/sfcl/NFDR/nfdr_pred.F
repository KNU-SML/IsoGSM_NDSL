#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NPES npes
#define LONF2 lonf2
#define LATG2 latg2
      subroutine nfdr_pred(iyear,month,iday,ihour,ahour,condir)
c
c  nfdr wrapper runs in mpi
c
      character*200 condir
c
#include <define.h>
#include <paramodel.h>
#ifndef RSM
#define LONF2 lonf2_
#define LATG2 latg2_
#define ILEN lonf_
#define JLEN latg_
#ifdef MP
#define MPGP2F mpgp2f
#define MPABORT mpabort
#include <commpi.h>
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else     /* MP */
#define MPABORT abort
#define LONF2S lonf2_
#define LATG2S latg2_
#endif      /* MP */
#else     /* RSM */
#define LONF2 igrd12_
#define LATG2 jgrd12_
#define ILEN igrd1_
#define JLEN jgrd1_
#ifdef MP
#define MPGP2F rmpgp2f
#define MPABORT rmpabort
#include <comrmpi.h>
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else     /* MP */ 
#define MPABORT abort 
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif      /* MP */
#endif      /* RSM */
c
#ifndef RSM
      real rproj,rorient,rtruth,rcotru,rdelx,rdely,
     1     rlat1,rlat2,rlon1,rlon2
      data rproj/-4./
      data rorient,rtruth,rcotru,rdelx,rdely,rlat1,rlat2,rlon1,rlon2
     *    /999.,999.,999.,999.,999.,999.,999.,0.,999./
#else
#include <rscomloc.h>
#endif
c
#include <comnfdr.h>
c#define DBGNFDR
c
      integer len
      parameter(len=ILEN*JLEN)
c
      real,allocatable :: dbgrd(:),rhgrd(:),pptgrd(:),cldgrd(:),
     *     wsgrd(:),pptdurgrd(:),tmaxgrd(:),tmingrd(:),
     *     rhmaxgrd(:),rhmingrd(:)
c
      allocate (dbgrd(len),rhgrd(len),pptgrd(len),cldgrd(len),
     *     wsgrd(len),pptdurgrd(len),tmaxgrd(len),tmingrd(len),
     *     rhmaxgrd(len),rhmingrd(len))
c
c  t2m2pm
c
#ifdef MP
      call MPGP2F(t2m2pm,LONF2S,LATG2S,dbgrd,LONF2,LATG2,1)
#else
      do n=1,len
        dbgrd(n)=t2m2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(dbgrd,1)
#else
      call rowsep(dbgrd)
#endif
c
c  rh2mtpm
c
#ifdef MP
      call MPGP2F(rh2m2pm,LONF2S,LATG2S,rhgrd,LONF2,LATG2,1)
#else
      do n=1,len
        rhgrd(n)=rh2m2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(rhgrd,1)
#else
      call rowsep(rhgrd)
#endif
c
c  precipitation
c
#ifdef MP
      call MPGP2F(ppt24,LONF2S,LATG2S,pptgrd,LONF2,LATG2,1)
#else
      do n=1,len
        pptgrd(n)=ppt24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(pptgrd,1)
#else
      call rowsep(pptgrd)
#endif
c
c  cloud cover
c
#ifdef MP
      call MPGP2F(cloud2pm,LONF2S,LATG2S,cldgrd,LONF2,LATG2,1)
#else
      do n=1,len
        cldgrd(n)=cloud2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(cldgrd,1)
#else
      call rowsep(cldgrd)
#endif
c
c  wind speed
c
#ifdef MP
      call MPGP2F(wind2pm,LONF2S,LATG2S,wsgrd,LONF2,LATG2,1)
#else
      do n=1,len
        wsgrd(n)=wind2pm(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(wsgrd,1)
#else
      call rowsep(wsgrd)
#endif
c
c  Precipitation duration
c
#ifdef MP
      call MPGP2F(prcpdur,LONF2S,LATG2S,pptdurgrd,
     1            LONF2,LATG2,1)
#else
      do n=1,len
        pptdurgrd(n)=prcpdur(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(pptdurgrd,1)
#else
      call rowsep(pptdurgrd)
#endif
c
c  max temperature
c
#ifdef MP
      call MPGP2F(tmax24,LONF2S,LATG2S,tmaxgrd,LONF2,LATG2,1)
#else
      do n=1,len
        tmaxgrd(n)=tmax24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(tmaxgrd,1)
#else
      call rowsep(tmaxgrd)
#endif
c
c  min temperature
c
#ifdef MP
      call MPGP2F(tmin24,LONF2S,LATG2S,tmingrd,LONF2,LATG2,1)
#else
      do n=1,len
        tmingrd(n)=tmin24(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(tmingrd,1)
#else
      call rowsep(tmingrd)
#endif
c
c  max relative humidity
c
#ifdef MP
      call MPGP2F(rhmax,LONF2S,LATG2S,rhmaxgrd,LONF2,LATG2,1)
#else
      do n=1,len
        rhmaxgrd(n)=rhmax(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(rhmaxgrd,1)
#else
      call rowsep(rhmaxgrd)
#endif
c
c  min relative humidity
c
#ifdef MP
      call MPGP2F(rhmin,LONF2S,LATG2S,rhmingrd,LONF2,LATG2,1)
#else
      do n=1,len
        rhmingrd(n)=rhmin(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(rhmingrd,1)
#else
      call rowsep(rhmingrd)
#endif
c
#ifdef MP
#ifdef RSM
      call rmpsynall
#else
      call mpsynall
#endif
#endif
c
      call nfdr(ILEN,JLEN,iyear,month,iday,ihour,ahour,
     *          rproj,rorient,rtruth,rcotru,
     *          rdelx,rdely,rlat1,rlat2,rlon1,rlon2,
     *          dbgrd,rhgrd,pptgrd,cldgrd,wsgrd,
     *          pptdurgrd,tmaxgrd,tmingrd,rhmaxgrd,rhmingrd,
     *          condir)
c
      deallocate (dbgrd,rhgrd,pptgrd,cldgrd,
     *     wsgrd,pptdurgrd,tmaxgrd,tmingrd,
     *     rhmaxgrd,rhmingrd)
c
#ifndef DBGNFDR
      ifhr=ahour-24.
      if(ifhr.lt.0) then
        print *,'ahour error in nfdr_pred'
        call abort
      endif
      ithr=ahour
      call nfdr_wriflx(88,iyear,month,iday,ihour,ifhr,ahour)
#ifdef MP
#ifdef RSM
      call rmpsynall
#else
      call mpsynall
#endif
#endif
c
#endif
      call nfdr_zer
c
      return
      end
