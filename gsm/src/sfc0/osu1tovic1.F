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
#define MTNVAR mtnvar
      subroutine osu1tovic1(sfcfcsin,ijdim,numsfcsin)
c
c fill vic1 type surface file records from osu1 type
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

#ifdef VICLSM1
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
          if(k.eq.1) then
             stc(ij,1,k) = tsea(ij,1)  ! the first node is ground surface
          else
             stc(ij,1,k)=sfcfcsin(ij,ind+k-1-1)
             stc(ij,1,k+1)=sfcfcsin(ij,ind+k-1)
          endif
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
#ifdef DBG
      open(166,file='CHECK_ini_temp.dat',status='unknown')
      do ij = 1, ijdim
         if(slmsk(ij,1).eq.1)then
            if(abs(stc(ij,1,1)-stc(ij,1,2)).gt.5.0) then
               write(166,1166) ij,stc(ij,1,1),stc(ij,1,2),stc(ij,1,3)
 1166          format('cell id=',i5,1x,'temp=',3f15.3)
               if(stc(ij,1,1).gt.stc(ij,1,2))then
                  stc(ij,1,2) = stc(ij,1,1) - 5.0
               else
                  stc(ij,1,2) = stc(ij,1,1) + 5.0
               endif
            endif
         endif
      enddo
      close(166)
#endif
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
c  fill 3-rd soil layer moisture and temperature
c
      do ij=1,ijdim
         do k=3,lsoil_
            smc(ij,1,k)=smc(ij,1,2)
         end do
         do k=4,nsoil_
            stc(ij,1,k)=stc(ij,1,3)
        enddo
      enddo
c
c  uustar, ffmm, ffhh, canopy are filled with reasonable constant
c
      alog30=log(30.)
      do ij=1,LONF2S*LATG2S
         uustar(ij,1)=1.
         ffmm(ij,1)=alog30
         ffhh(ij,1)=alog30
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
c  fill prcp by zero
c
      do ij=1,LONF2S*LATG2S
         prcp  (ij,1)=0.
      enddo
c
c  fill vfrac vegtyp, albdeo, albedo fraction with 1.e30.  
c       vroot, binf, Ds, Dsm, Ws,cef, expt, kest, dph, bub, qrt, 
c       bkd, sld, wcr, wpw, smr, smx, dphn, smxn, expn, bubn, 
c       alpn, betn, gamn, flai, silz, snwz, sic, csno, rsno,
c       tsf, tpk, sfw, pkw, lstsn
c
c  These fields will be replaced with climatology field by sfcmrg call
c
      undef=1.e30
      do ij=1,LONF2S*LATG2S
         vfrac(ij,1)=undef
         vtype(ij,1)=undef
         do l=1,2
            facalf(ij,1,l)=undef
         enddo
         do l=1,lalbd_
            albedo(ij,1,l)=undef
         enddo
c
         binf(ij,1) = undef
         Ds  (ij,1) = undef
         Dsm (ij,1) = undef
         Ws  (ij,1) = undef
         cef (ij,1) = undef
         silz(ij,1) = undef
         snwz(ij,1) = undef
         csno(ij,1) = undef
         rsno(ij,1) = undef
         tsf (ij,1) = undef
         tpk (ij,1) = undef
         sfw (ij,1) = undef
         pkw (ij,1) = undef
         lstsn(ij,1)= undef
         do k = 1, lsoil_
            vroot(ij,1,k)= undef
            expt(ij,1,k) = undef
            kest(ij,1,k) = undef
            dph (ij,1,k) = undef
            bub (ij,1,k) = undef
            qrt (ij,1,k) = undef
            bkd (ij,1,k) = undef
            sld (ij,1,k) = undef
            wcr (ij,1,k) = undef
            wpw (ij,1,k) = undef
            smr (ij,1,k) = undef
            smx (ij,1,k) = undef
            sic (ij,1,k) = undef
         end do
         do k = 1, nsoil_
            dphn(ij,1,k) = undef
            smxn(ij,1,k) = undef
            expn(ij,1,k) = undef
            bubn(ij,1,k) = undef
            alpn(ij,1,k) = undef
            betn(ij,1,k) = undef
            gamn(ij,1,k) = undef
         end do
c
      enddo
c
#endif
      return
      end
