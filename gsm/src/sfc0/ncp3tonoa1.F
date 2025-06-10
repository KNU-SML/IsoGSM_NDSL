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
c
c
      subroutine ncp3tonoa1(sfcfcsin,ijdim,numsfcsin)
c
c
c                                                    (see sfcfld.F           )
c fill noa1 type surface file records from ncp3 type (200509T574L64 1760x880)
c
c
      implicit none
c
#include <paramodel.h>
#include <comsfc.h>
c
      integer ijdim,numsfcsin
      real sfcfcsin(ijdim,numsfcsin)
c
      integer ind,ij,k,l
      real alog30,undef
c
c
c
#ifdef NOALSM1
c
c tsea
c
      ind=3
      do ij=1,ijdim
        tsea(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c smc
c
      ind=37
      do k=1,4
        do ij=1,ijdim
          smc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
c
c snow
c
      ind=4
      do ij=1,ijdim
        sheleg(ij,1)=max(sfcfcsin(ij,ind)*100,0.)
      enddo
c
c stc
c
      ind=33
      do k=1,4
        do ij=1,ijdim
          stc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
c
c tg3
c
      ind=5
      do ij=1,ijdim
        tg3(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c zorl
c
      ind=6
      do ij=1,ijdim
        zorl(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c cv
c
      undef=0.0
      do ij=1,ijdim
        cv(ij,1)=undef
      enddo
c
c cvb
c
      do ij=1,ijdim
        cvb(ij,1)=undef
      enddo
c
c cvt
c
      do ij=1,ijdim
        cvt(ij,1)=undef
      enddo
c
c albedo
c
      ind=7
      do k=1,4   ! ? missing until Aug 2006
         do ij=1,ijdim
           albedo(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      enddo      ! ? missing until Aug 2006
c
c slmsk
c
      ind=1
      do ij=1,ijdim
        slmsk(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c vegetation cover
c
      ind=11
      do ij=1,ijdim
        vfrac(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c canopy
c
      ind=12
      do ij=1,ijdim
        canopy(ij,1)=max(sfcfcsin(ij,ind),0.)
      enddo
c
c f10m
c
      ind=13
      do ij=1,ijdim
        f10m(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c vegitation type
c
      ind=16
      do ij=1,ijdim
        vtype(ij,1)=max(sfcfcsin(ij,ind),0.)
      enddo
c
c soil type
c
      ind=17
      do ij=1,ijdim
        stype(ij,1)=max(sfcfcsin(ij,ind),0.)
      enddo
c
c albedo fraction type
c
      ind=18
      do k=1,2
        do ij=1,ijdim
          facalf(ij,1,k)=max(sfcfcsin(ij,ind+k-1),0.)
        enddo
      enddo
c
c ustar
c
      ind=20
      do ij=1,ijdim
        uustar(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c ffmm
c
      ind=21
      do ij=1,ijdim
        ffmm(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c ffhh
c
      ind=22
      do ij=1,ijdim
        ffhh(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c prcp
c
      ind=26
      do ij=1,ijdim
        prcp(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c srflag
c
      ind=27
      do ij=1,ijdim
        srflag(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c snwdph
c
      ind=28
      do ij=1,ijdim
        snwdph(ij,1)=max(sfcfcsin(ij,ind)*100,0.)
      enddo
c
c slc
c
      ind=41
      do k=1,4
        do ij=1,ijdim
          slc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
c
c shdmin
c
      ind=29
      do ij=1,ijdim
        shdmin(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c shdmax
c
      ind=30
      do ij=1,ijdim
        shdmax(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c slope
c
      ind=31
      do ij=1,ijdim
        slope(ij,1)=sfcfcsin(ij,ind)
      enddo
c
c snoalb
c
      ind=32
      do ij=1,ijdim
        snoalb(ij,1)=max(sfcfcsin(ij,ind),0.)
      enddo
c
#endif
      return
      end
