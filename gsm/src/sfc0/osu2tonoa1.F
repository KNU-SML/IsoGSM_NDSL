#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define NSOIL nsoil
      subroutine osu2tonoa1(sfcfcsin,ijdim,numsfcsin)
c
c fill noa1 type surface file records from osu1 type
c
#include <paramodel.h>
#include <comsfc.h>
c
      integer ijdim,numsfcsin
      real sfcfcsin(ijdim,numsfcsin)
c
      integer ind,ij,k,l
      real alog30,undef
#ifdef NOALSM1
c
c tsea
c
      ind=1
      do ij=1,ijdim
        tsea(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c smc
c
      do k=1,2
        do ij=1,ijdim
          smc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
      ind=ind+2
c
c snow
c
      do ij=1,ijdim
        sheleg(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c stc
c
      do k=1,2
        do ij=1,ijdim
          stc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
      ind=ind+2
c
c tg3
c
      do ij=1,ijdim
        tg3(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c zorl
c
      do ij=1,ijdim
        zorl(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c cv
c
      do ij=1,ijdim
        cv(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c cvb
c
      do ij=1,ijdim
        cvb(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c cvt
c
      do ij=1,ijdim
        cvt(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c albedo
c
      do k=1,lalbd_
        do ij=1,ijdim
          albedo(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
        ind=ind+1
      enddo
c
c slmsk
c
      do ij=1,ijdim
        slmsk(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c vegetation cover
c
      do ij=1,ijdim
        vfrac(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c canopy
c
      do ij=1,ijdim
        canopy(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c f10m
c
      do ij=1,ijdim
        f10m(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c vegitation type
c
      do ij=1,ijdim
        vtype(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c soil type
c
      do ij=1,ijdim
        stype(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c albedo fraction type
c
      do k=1,2
        do ij=1,ijdim
          facalf(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
        ind=ind+1
      enddo
c
c ustar
c
      do ij=1,ijdim
        uustar(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c ffmm
c
      do ij=1,ijdim
        ffmm(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c ffhh
c
      do ij=1,ijdim
        ffhh(ij,1)=sfcfcsin(ij,ind)
      enddo
c
      if(ind.ne.numsfcsin) then
         print *,'counting error in osu2tonoa1'
         call abort
      endif
c
c  fill 3-rd and 4-th soil layer moisture and temperature
c
      do ij=1,LONF2S*LATG2S
        do k=3,lsoil_
          smc(ij,1,k)=smc(ij,1,2)
          stc(ij,1,k)=stc(ij,1,2)
        enddo
      enddo
c
c  initialize snwdph from sheleg
c
      do ij=1,ijdim
        snwdph(ij,1)=sheleg(ij,1)*5.0
      enddo
c
c  initialize srflag from tsea (t850 is not available)
c
      print *,'initialize srflag-use tsea as surrogate for t850'
      do ij=1,ijdim
        srflag(ij,1)=0.
        if(tsea(ij,1).le.273.16) srflag(ij,1)=1.
      enddo
c
c  initialize lsoil layer slc from smc and stc
c
      call getslc(slmsk,stype,stc,smc,ijdim,lsoil_,slc)
c
c  fill prcp by zero
c
      do ij=1,LONF2S*LATG2S
        prcp  (ij,1)=0.
      enddo
c
c  fill vfrac vegtyp, stype, albdeo, albedo fraction with 1.e30.  
c  These fields will be replaced with climatology field in sfc program.
c
      undef=1.e30
      do ij=1,LONF2S*LATG2S
        shdmin(ij,1)=undef
        shdmax(ij,1)=undef
        snoalb(ij,1)=undef
        slope (ij,1)=undef
      enddo
c
#endif
      return
      end
