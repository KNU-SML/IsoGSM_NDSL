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
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NPES npes
#define NTOTAL ntotal
#define NTRAC ntrac
#define PI pi
#include <define.h>
      subroutine grrad1(lons2,
#ifdef SCM
     &                  dgrs,
#else
     &                  dgrs,plamgr,pphigr,
#endif
     &                  ugrs,vgrs,tgrs,qgrs,pgr,
     1                  albedr,slmskr,rlon,rlat,
     2                  tsear,shelgr,tgr,cvr,cvtr,cvbr,rhcl,
#ifndef SW99UPDATE
     3                  ozonea,albdoa,cldary,cldtot,cldcnv,
#else
     3                  ozonea,albdoa,cldary,cldtot,cldcnv,rhrh,
#endif
     4                  cldsa,mtopa,mbota,
     5                  rrs2,lat,latco,latrue,istrat,
     6                  kalb,jo3,slag,rsin1,rcos1,rcos2,
     7                  fjd,dlt,jsno,workr,lworkr
#ifdef DG
     8                  ,cldt,cldc
#endif
#ifdef DGZ
     9                  ,sfcp
#endif
#ifdef CONDENS
     +                  ,qcicnv, qrscnv
     +                  ,qcilrg, qrslrg
#endif
#ifdef SW99UPDATE
     b                  ,icwp,qcis
#endif
     b                 )
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
c
c     ************************************************************
c     *  added accumulation of clds and convective cloud in dg3  *
c     *                                        k.a.c sept 1994   *
c     *  changed  h,m,l calculation in cldjms (removed facv)     *
c     *     and added proper total cloud calculation             *
c     *  changed avecld calc in cldiag (used cldary) and         *
c     *     used total cloud calculated in cldjms                *
c     *                                         k.a.c. nov94     *
c     *  interpolate o3 profile to each gridpoint, ie use        *
c     *     proper surface pressure                              *
c     *                                         k.a.c. dec94     *
c     *  fix pl1 for operations, where dgz is on and dg3 is off, *
c     *     ....note dg is on if either dgz or dg3 is on         *
c     *                                         k.a.c. jan94     *
c     ************************************************************
c
c        updates made to add oceanic stratus and to fix conv cloud..
c                     to gloopr - ivv(2),ibl are set=1....
c                     to gloopr - set min q to 1.e-10,rather than 1.e-6
c                                 to anticipate avoiding cld creation
c                                 in extremely dry,cold (winter) regions
c                                 where 1.e-6 could imply hi valu of rh
c                     to cldjms - multitude of changes
c        updates made to add grid point diagnostics ..k.a.c...sep 91
c                     to gloopr -
c        updates made to fix sw approx              ..k.a.c...nov 91
c                     to coszmn
c        updates made to pass and receive sib data  ..k.a.c...mar 92
c                     to gloopr -
c        updates made to fix sw rad diagnostics     ..k.a.c...jun 92
c                     proper diurnal weighting
c                     to gloopr and coszmn
c        updates made to calculate clear-sky "on-the-fly" kac aug 92
c                     to gloopr,radfs,fst,spa,lwr,swr
c                     ...for cloud forcing....
c        updates made for the completely new cloud routine (cldjms),use
c                     flag ivva to control vertical velocity adj.
c                     for low cld (=0: without,  =1: with)
c                     use flag invr to control lapse rate inversion
c                     type of cld (=0: without,  =1: with)
c                     to gloopr and radfs ...y.h.           ...dec92
c        updates made to call cld optical property routine (cldprp),
c                     to give cld emissivity, optical depth, layer
c                     reflectance and transmitance
c                     to gloopr and radfs ...y.h.           ...feb93
c
c       updates made to allow tuned cld-rh data to be used..ctune
c                   to gloopr and cldjms          ..k.a.c...dec 92
c                   spatial interpolation of tables ........may93
c                   use only 1 set of tuning tables for all fcst hrs,
c                    the tuning of the 24hr fcst       .....jan94
c                    old code used 6 tables..see ckc94 .....feb94
c                   since tuning done for h,m,l cld, vertically
c                    blend the relations at old hml bdries..jan94
c        updates made to change definition of h,m,l domains..
c          to gloopr,gcljms,cldjms, cldprp   k.a.c...dec92 + aug93
c          wmo aerosols distributions, and b.p.breiglebs surface
c          albedo scheme.            ..............y.h...sep95
c             kalb is the control flag for surface albedo
c             kalb=0 use the old matthews data, =1 use the new scheme.
c             jo3 is the control flag for ozone climatology
c             jo3=0 use the old gfdl data, =1 use the new nasa data
#ifdef SW99UPDATE
c        recoded the default part of cloud scheme patterned after
c          ccm3 (kiehl et al. 1998,j.clim; and 1994,jgr)... y.h. aug98
c        updates made to opac aerosol algorithm (1998)..... y.h. sep99
c        updates made to cloud overlapping and moorthi's liquid
c          water cloud fraction scheme ............y.h...oct00
#endif
c
#ifndef SW99UPDATE
#ifndef ICECLOUD
#include <rdparm.h>
#else
#include <rdparm8.h>
#endif
#else
#include <rdparm99.h>
#endif
c
#ifndef RSM
#include <comfver.h>
#define LONLENS lons2
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#else
#include <rscomver.h>
#ifdef MP
#include <comrmpi.h>
#ifdef RSMVECTORIZE
#define LONF2S igrd12p_*jgrd12p_
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONF2S igrd12p_
#define LONLENS lonlen(mype)*2
#endif
#define LATG2S  jgrd12p_
#else
#undef RSMVECTORIZE
#define LONF2S igrd12_
#define LONLENS igrd12_
#define LATG2S  jgrd12_
#endif
#endif
c
      parameter (mcld=3,nseal=2,nbin=100,nlon=2,nlat=4)
c
      dimension dgrs(LONF2S,levs_)
#ifndef SCM
      dimension plamgr(LONF2S),pphigr(LONF2S)
#endif
      dimension ugrs(LONF2S,levs_),
     *          vgrs(LONF2S,levs_)
      dimension tgrs(LONF2S,levs_),
     *          qgrs(LONF2S,levh_)
      dimension pgr(LONF2S)
      dimension workr(LONF2S,lworkr)
#ifndef RSM
      dimension rrs2(LATG2S)
#endif
      dimension rlon(LONF2S),rlat(LONF2S)
      dimension albedr(LONF2S),slmskr(LONF2S)
      dimension vvel(LONF2S,levs_),
     *          ssnow(LONF2S)
c
c   cldary contains multi layers of cloud
c
#ifndef SW99UPDATE
      dimension cldary(LONF2S,levs_)
#else
      dimension cldary(LONF2S,levs_),
     *          rhrh (LONF2S,levs_)
#endif
      dimension cldcnv(LONF2S,levs_),
     *          cldtot(LONF2S,levs_)
#ifdef DGZ
      dimension sfcp(LONF2S)
#endif
#ifdef DG
      dimension icnv(LONF2S)
      dimension cldc(LONF2S,levs_),
     *          cldt(LONF2S,levs_)
#endif
#ifdef CONDENS
      dimension qcicnv(LONF2S,levs_),qrscnv(LONF2S,levs_)
      dimension qcilrg(LONF2S,levs_),qrslrg(LONF2S,levs_)
#endif
#ifdef SW99UPDATE
      dimension qcis(LONF2S,levs_)
#endif
      dimension ozonea(LONF2S,levs_),
     *          albdoa(LONF2S),
     *          cldsa(LONF2S,4),
     *          mtopa(LONF2S,3),
     *          mbota(LONF2S,3),
     *          wrkema(LONF2S),
     *          tsear(LONF2S),shelgr(LONF2S),
     *          cvr(LONF2S),cvtr(LONF2S),
     *          cvbr(LONF2S),tgr(LONF2S)
chybrid
      dimension sihyb(LONF2S,levp1_),cihyb(LONF2S,levp1_)
      dimension slhyb(LONF2S,levs_),clhyb(LONF2S,levs_)
      dimension delhyb(LONF2S,levs_),slkhyb(LONF2S,levs_)
c
      save invr,ivva,rhmax,xlabdy,xlobdy,xlim,ko3cli
c
c  array added for rh-cl calculation
c  indices for lon,lat,cld type(l,m,h), land/sea respectively
c  nlon=1-2, for eastern and western hemispheres
c  nlat=1-4, for 60n-30n,30n-equ,equ-30s,30s-60s
c  land/sea=1-2 for land(and seaice),sea
c
      dimension rhcl (nbin,nlon,nlat,mcld,nseal)
      dimension rhcla(nbin,nlon,mcld)
      dimension rhcld(LONF2S,nbin,mcld)
      dimension xlabdy(3),xlobdy(3)
c
c  xlabdy = lat bndry between tuning regions,+/- xlim for transition
c  xlobdy = lon bndry between tuning regions
c
      data xlabdy / 30.e0 , 0.e0 , -30.e0 /
      data xlobdy / 0.e0 , 180.e0 , 360.e0 /
      data xlim / 5.e0 /
c
c  initial rh crit. set 1 for ocean, set 2 for land.
c  invr=0 no lapse rate inversion type cld, =1 wiht it
c  ivva=0 no vertical velocity adj. for low cld, =1 with adj.
c
      data rhmax/1.00e0/, invr/1/, ivva/1/
c
c   ko3cli=0,1 for gfdl,nasa climo
c     only for ko3=1(meaning o3 is input to radfs)..see grrad2
c
      ko3cli=jo3
#ifndef RSM
#ifdef DGZ
      do 250 j=1,LONLENS
        sfcp(j) = pgr(j)
  250 continue
#endif
#else
      do 250 j=1,LONLENS
        pgr(j) = exp(pgr(j))
  250 continue
#endif
#ifdef DG
      do 260 j=1,LONLENS
        icnv(j) = 0.e0
  260 continue
#endif
#ifdef HYBRID
      do j=1,lons2
        do k=1,levp1_
          sihyb(j,k)=ak5(levp1_-k+1)/pgr(j)+bk5(levp1_-k+1)
          cihyb(j,k)=1-sihyb(j,k)
        enddo
        do k=1,levs_
          slhyb(j,k)=0.5*(sihyb(j,k)+sihyb(j,k+1))
          delhyb(j,k)=sihyb(j,k)-sihyb(j,k+1)
          clhyb(j,k)=1-slhyb(j,k)
        enddo
      enddo
#else
      do j=1,lons2
        do k=1,levp1_
          sihyb(j,k)=si(k)
          cihyb(j,k)=ci(k)
        enddo
        do k=1,levs_
          slhyb(j,k)=sl(k)
          delhyb(j,k)=del(k)
          clhyb(j,k)=cl(k)
        enddo
      enddo
#endif
c
c   generalized the computation and suitable for all grid ---
c                                             by h.-m.juang
      do 290 i=1,LONLENS
c
      isla = 1
      if (slmskr(i).lt.1.e0) isla = 2
      xlatpt = rlat(i) * 180.e0 / pi_
c
c  get rh-cld relation for this lat
c
      kla = 4
      do k=1,3
       if (xlatpt.gt.xlabdy(k)) then
        kla = k
        go to 215
       end if
      enddo
  215 continue
      klap=0
      do k=1,3
       xlnn = xlabdy(k)+xlim
       xlss = xlabdy(k)-xlim
       if (xlatpt.lt.xlnn.and.xlatpt.gt.xlss) then
         kla  = k
         klap = k+1
         go to 217
       endif
      enddo
  217 continue
      if( klap .eq. 0 ) then
       do kc=1,mcld
        do lo=1,nlon
         do nbi=1,nbin
            rhcla(nbi,lo,kc) = rhcl(nbi,lo,kla,kc,isla)
         enddo
        enddo
       enddo
      else
c
c  linear transition between latitudinal regions...
c
       do kc=1,mcld
        do lo=1,nlon
         do nbi=1,nbin
            rhcla(nbi,lo,kc) =
     1         (rhcl(nbi,lo,kla,kc,isla)-rhcl(nbi,lo,klap,kc,isla))
     2       * (xlatpt-xlss)/(xlnn-xlss) + rhcl(nbi,lo,klap,kc,isla)
          enddo
        enddo
       enddo
      endif
c
c  get rh-cld relation for this lon
c
      xlonpt = rlon(i) * 180.e0 / pi_
       lo=1
       if( xlonpt .gt. 180.0 ) lo=2
       ikn = 0
        do k=1,3
         diflo = abs(xlonpt-xlobdy(k))
         if (diflo.lt.xlim) then
          ikn = k
          ilft = lo
          irgt = ilft + 1
          if (irgt.gt.nlon) irgt = 1
          xlft = xlobdy(ikn) - xlim
          xrgt = xlobdy(ikn) + xlim
          go to 219
         endif
        enddo
 219   continue
       if( ikn.eq.0 ) then
         do k=1,mcld
          do nbi=1,nbin
            rhcld(i,nbi,k) = rhcla(nbi,lo,k)
          enddo
         enddo
       else
         do k=1,mcld
          do nbi=1,nbin
            rhcld(i,nbi,k) =
     1     (rhcla(nbi,ilft,k)-rhcla(nbi,irgt,k))
     2      * (xlonpt-xrgt)/(xlft-xrgt)+rhcla(nbi,irgt,k)
          enddo
         enddo
       endif
c
  290   continue
#ifndef SCM
c
c   get mean zenith angle for this dtswav-both nh and sh
c   start radfs set-up for both hemispheres
c   get vertical motion (cb/sec) in vvel
c
#ifdef RSMVECTORIZE
#define OMEGAS omegasv
#else
#define OMEGAS omegas
#endif
      call OMEGAS(LONLENS,LONF2S,levs_,
     &            pphigr(1),plamgr(1),workr,ugrs(1,1),vgrs(1,1),
#ifndef RSM
     2            dgrs(1,1),delhyb,rrs2(lat),vvel,pgr(1),slhyb)
#else
     2            dgrs(1,1),delhyb,rrs2     ,vvel,pgr(1),slhyb)
#endif
#endif
c
c  get model diagnosed clds
c
#ifndef SW99UPDATE
#ifndef CLDSLINGO
#define CLDJMS cldjms
#else
#define CLDJMS cldjmsn
#endif
      call CLDJMS(lons2,LONF2S,levs_,nbin,mcld,
     1            pgr(1),qgrs(1,1),tgrs(1,1),vvel,
     2            cvr(1),cvtr(1),cvbr(1),
     3            sihyb,slhyb,    slmskr(1),
     4            cldsa(1,1),mtopa(1,1),mbota(1,1),
     5            cldary(1,1),ivva,invr,rhmax,
#ifdef DBG
     6            rlat(1),rhcld,istrat,ntotal_,lat)
#else
     6            rlat(1),rhcld,istrat,ntotal_)
#endif
c
cdg3     unpack cldamt and conv cldamt from cldary..(stratus is +2)
cdg3      radiation sees stratiform or convective...not merged.....
cdg3      so cldtot reflects this below
cdg3  note:cldary 2-4 digits to left of decimal=cv cloud
cdg3       cldary 1 digit to left of decimal+fractional part=strat cld
cdg3        get stratiform cld into cldtot, conv cld in cldcnv
cdg3        ..anvil ci is stratiform, so place the 1 lyr in cnv also
cdg3        ..icnv = cv cloud top layer
c
       do 251 k=1,levs_
        do 251 i=1,LONLENS
         cldtot(i,k) = amod (cldary(i,k),2.e0)
         cldcnv(i,k) = float(int(cldary(i,k))/10)*1.e-3
  251  continue
#else
      call cldjms00(lons2,LONF2S,levs_,nbin,mcld,
     1            pgr(1),qgrs(1,1),tgrs(1,1),vvel,
     2            cvr(1),cvtr(1),cvbr(1),
     3            sihyb,slhyb,    slmskr(1),
     4            cldsa(1,1),mtopa(1,1),mbota(1,1),
     5            cldtot(1,1),cldcnv(1,1),ivva,invr,icwp,qcis,rhrh,
#ifdef DBG
     6            rlat(1),rhcld,istrat,ntotal_,lat)
#else
     6            rlat(1),rhcld,istrat,ntotal_)
#endif
#endif
c
cdg3   anvil ci is stratiform, so place the 1 lyr in cnv also
cdg3   icnv = cv cloud top layer
c
#ifdef DG
       do k=1,levs_
        do i=1,LONLENS
         cldt(i,k) = cldtot(i,k)
         cldc(i,k) = cldcnv(i,k)
         if(cldc(i,k).gt.0.e0) then
          cldt(i,k) = cldc(i,k)
          icnv(i) = k
         end if
        end do
       end do
c
cdg3   anvil ci placed into layer above conv cld top
c
       do 252 i=1,LONLENS
        k = icnv(i)
        if (k.gt.1.and.k.lt.levs_) then
         cldc(i,k+1) = cldt(i,k+1)
        end if
  252  continue
#endif
c
c  compute mean cloud diagnostics + h,m,l,total cloud
c     call cvdiag(lons2,avecv(1,lat  ),avecv(1,latco),
c    1                         cvr(1),cvtr(1),cvbr(1))
c     call cldiag(lons2,avecld(1,lat  ),cldl(1,lat  ),
c    1                  avecld(1,latco),cldl(1,latco),
c    2                         cldsa(1,1),cldary(1,1))
c
#ifndef O3C
      if(thour.eq.0.or.ntrac_.le.1) then
#endif
        if (ko3cli.eq.0)
     *    call ozon2d(LONLENS,LONF2S,levs_,
     *                slhyb,ozonea(1,1),
     *                pgr(1),rlat(1),rsin1,rcos1,rcos2)
        if (ko3cli.eq.1)
     *    call oz2d(LONLENS,LONF2S,levs_,
     *              slhyb,ozonea(1,1),
     *              pgr(1),rlat(1))
#ifndef O3C
      endif
#endif
c
      do 340 i=1,LONLENS
        if(slmskr(i).eq.2.0e0) then
           tsear(i) = min(tsear(i),271.2e0)
        else if(slmskr(i).eq.0.0e0) then
           tsear(i) = max(tsear(i),271.21e0)
        else if(slmskr(i).eq.1.0e0 .and. shelgr(i).gt.0.0e0)
     1  then
           tsear(i) = min(tsear(i),273.16e0)
        endif
  340 continue
      if (kalb .ge. 1) then
        do i=1,LONLENS
          albdoa(i)=0.
        enddo
        return
      endif
c
c   the following determines surface albedo (albdoa),where snow exists.
c
      do 380 i=1,LONLENS
        wrkema(i) = 0.e0
        ssnow(i) = shelgr(i) * 0.1e0
  380 continue
#ifndef RSM
      call albsno(LONLENS,latrue,jsno,
#else
      call albsno(LONLENS,lat,jsno,
#endif
     1            albdoa(1),rlat(1),albedr(1),slmskr(1),
     2            ssnow, tgr(1), tgrs(1,1) )
c
      return
      end
