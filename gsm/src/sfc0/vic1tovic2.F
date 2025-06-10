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
#define NSOIL nsoil
#define MSUB  msub
#define MTNVAR mtnvar
c
      subroutine vic1tovic2(sfcfcsin,ijdim,numsfcsin)
c
c fill vic2 type surface file records from noa1 type
c
#include <paramodel.h>
#include <comsfc.h>
c
      integer ijdim,numsfcsin
      real sfcfcsin(ijdim,numsfcsin)
c
      integer ind,ij,k,l,nv
      real undef
c
#ifdef VICLSM2
c
      integer kslmb_, nslmb_
      parameter (kslmb_=lsoil_*msub_, nslmb_=nsoil_*msub_)
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
      do k=1, lsoil_
        do ij=1,ijdim
          smc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
      ind=ind+lsoil_
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
      do k=1,nsoil_
        do ij=1,ijdim
          stc(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
      ind=ind+nsoil_
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
c albedo    !vic1 kalbd=4
c
      do k = 1, 4
         do ij=1,ijdim
            albedo(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+4
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
c soil type (in VIC no soil type)
c
      ind=ind+1
c
c albedo fraction type
c
      do k=1,2
        do ij=1,ijdim
          facalf(ij,1,k)=sfcfcsin(ij,ind+k-1)
        enddo
      enddo
      ind=ind+2
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
      ind=ind+1
c
c  prcp
c
      do ij=1,ijdim
        prcp(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  srflag
c
      do ij=1,ijdim
        srflag(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  binf
c
      do in=1,ijdim
         binf(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  Ds
c
      do in=1,ijdim
         Ds(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  Dsm
c
      do in=1,ijdim
         Dsm(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  Ws
c
      do in=1,ijdim
         Ws(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  cef
c
      do in=1,ijdim
         cef(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind=ind+1
c
c  expt
c
      do k = 1, lsoil_
         do in=1,ijdim
            expt(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  kst
c
      do k = 1, lsoil_
         do in=1,ijdim
            kst(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  dph
c
      do k = 1, lsoil_
         do in=1,ijdim
            dph(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  bub
c
      do k = 1, lsoil_
         do in=1,ijdim
            bub(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  qrt
c
      do k = 1, lsoil_
         do in=1,ijdim
            qrt(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  bkd
c
      do k = 1, lsoil_
         do in=1,ijdim
            bkd(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  sld
c
      do k = 1, lsoil_
         do in=1,ijdim
            sld(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  wcr
c
      do k = 1, lsoil_
         do in=1,ijdim
            wcr(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  wpw
c
      do k = 1, lsoil_
         do in=1,ijdim
            wpw(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  smr
c
      do k = 1, lsoil_
         do in=1,ijdim
            smr(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  smx
c
      do k = 1, lsoil_
         do in=1,ijdim
            smx(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  dphn
c
      do k = 1, nsoil_
         do in=1,ijdim
            dphn(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+nsoil_
c
c  smxn
c
      do k = 1, nsoil_
         do in=1,ijdim
            smxn(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+nsoil_
c
c  expn
c
      do k = 1, nsoil_
         do in=1,ijdim
            expn(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+nsoil_
c
c  bubn
c
      do k = 1, nsoil_
         do in=1,ijdim
            bubn(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+nsoil_
c
c  alpn
c
      do k = 1, nsoil_
         do in=1,ijdim
            alpn(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+nsoil_
c
c  flai
c
      ind = ind + 1
c
c  silz
c
      do in=1,ijdim
         silz(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  snwz
c
      do in=1,ijdim
         snwz(ij,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  sic
c
      do k = 1, lsoil_
         do in=1,ijdim
            msic(ij,1,k)=sfcfcsin(ij,ind+k-1)
         enddo
      end do
      ind=ind+lsoil_
c
c  csno
c
      do in=1,ijdim
         csno(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  rsno
c
      do in=1,ijdim
         rsno(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  tsf
c
      do in=1,ijdim
         tsf(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  tpk
c
      do in=1,ijdim
         tpk(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  sfw
c
      do in=1,ijdim
         sfw(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  pkw
c
      do in=1,ijdim
         pkw(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1
c
c  lstsn
c
      do in=1,ijdim
         lstsn(ij,1,1)=sfcfcsin(ij,ind)
      enddo
      ind = ind + 1

      if(ind.ne.numsfcsin) then
         print *,'counting error in vic1tovic2'
         call abort
      endif
c
c  fill subgrid variables
c
      do ij = 1, ijdim
         if(msub_.gt.1) then
            do nv = 1, msub_
               mcnp(ij,1,nv) = canopy(ij,1)
               msno(ij,1,nv) = msno(ij,1,1)
               csno(ij,1,nv) = csno(ij,1,1)
               rsno(ij,1,nv) = rsno(ij,1,1)
               tsf(ij,1,nv)  = tsf(ij,1,1)
               tpk(ij,1,nv)  = tpk(ij,1,1)
               sfw(ij,1,nv)  = sfw(ij,1,1)
               pkw(ij,1,nv)  = pkw(ij,1,1)
               lstsn(ij,1,nv)= lstsn(ij,1,1)
            end do
         endif
      enddo
c
c  fill vroot, binf, nveg,flai,msno,msmc,msic,mstc
c
c  These fields will be replaced with climatology field in sfc program.
c
      undef=1.e30
      do ij=1,LONF2S*LATG2S
         do nv = 1, msub_
            mvfr(ij,1,nv)=undef
            mvty(ij,1,nv)=undef
         enddo
c
         nveg(ij,1) = undef
         do k = 1, kslmb_
            vroot(ij,1,k) = undef
            msmc (ij,1,k) = undef
            msic (ij,1,k) = undef
         enddo
         do k = 1, nslmb_
            mstc(ij,1,k) = undef
         enddo
c
      enddo
c
#endif
      return
      end
