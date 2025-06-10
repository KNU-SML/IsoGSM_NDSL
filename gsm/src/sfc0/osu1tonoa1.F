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
      subroutine osu1tonoa1(sfcfcsin,ijdim,numsfcsin)
c
c fill noa1 type surface file records from osu1 type
c
      implicit none
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
      do ij=1,ijdim
        albedo(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c slmsk
c
      do ij=1,ijdim
        slmsk(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c plantr to be discarded
c
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
c
      if(ind.ne.numsfcsin) then
         print *,'counting error in osu1toosu2'
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
c  default soil type assumed as sandy loam (coarse-medium)
c
      do ij=1,ijdim
        stype(ij,1)=4
      enddo
      call getslc(slmsk,stype,stc,smc,ijdim,lsoil_,slc)
c
c  uustar, ffmm, ffhh, canopy are filled with reasonable constant
c
      alog30=log(30.)
      do ij=1,LONF2S*LATG2S
        uustar(ij,1)=1.
        ffmm(ij,1)=alog30
        ffhh(ij,1)=alog30
        prcp  (ij,1)=0.
      enddo
c
c  fill vfrac vegtyp, stype, albdeo, albedo fraction with 1.e30.  
c  These fields will be replaced with climatology field by sfcmrg call
c
      undef=1.e30
      do ij=1,LONF2S*LATG2S
        vfrac(ij,1)=undef
        vtype(ij,1)=undef
        stype(ij,1)=undef
        do l=1,2
          facalf(ij,1,l)=undef
        enddo
        do l=1,lalbd_
          albedo(ij,1,l)=undef
        enddo
        shdmin(ij,1)=undef
        shdmax(ij,1)=undef
        snoalb(ij,1)=undef
        slope (ij,1)=undef
      enddo
c
#endif
      return
      end
