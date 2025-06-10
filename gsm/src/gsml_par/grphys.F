#define BGF bgf
#define BORDER border
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
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define PI pi
#define QMIN qmin
#include <define.h>
c
c  rsm variables:
c
c    ddgrd = xxgrs
c    pxgrd = pllamgr
c    pygrd = plphigr
c    uugrd = uugrs
c    vvgrd = vvgrs
c    ttgrd = tvgrs
c    qqgrd = rqgrs
c    ppgrd = plgr 
c    
      subroutine grphys(lons2,lats2,
#ifdef SCM
     &                  uugrs,vvgrs,plgr,tvgrs,rqgrs,xxgrs,
#else
     &                  pllamgr,plphigr,uugrs,vvgrs,plgr,
     &                  tvgrs,rqgrs,xxgrs,
#endif
#ifdef DG3
     &                  gda,
#endif
#ifdef DG
     &                  cldt,clcv,
#endif
#ifdef DGZ
     &                  sfcp,
#endif
#ifdef CONDENS
     +                  taucld,cldwp,cldip,
#endif
     &                  lat,latrue)
#ifdef CRAY_THREAD
cfpp$ expand(znlacm,accdia)
#endif
c     ************************************************************
c     *  added accumulation of clds and convective cloud in dg3  *
c     *                                        k.a.c sept 1994   *
c     *  f3d added for clouds..mi new code=f94/source2/diagnew   *
c     *                               b katz + k.a.c oct  1994   *
c     *  changed  h,m,l calculation in cldjms (removed facv)     *
c     *     and added proper total cloud calculation             *
c     *  changed avecld calc in cldiag (used cldary) and         *
c     *     used total cloud calculated in cldjms                *
c     *  changed kenpts to store total cloud and all lyrs        *
c     *     of cloud                                             *
c     *                                         k.a.c. nov94     *
c     *  interpolate o3 profile to each gridpoint, ie use        *
c     *     proper surface pressure                              *
c     *                                         k.a.c. dec94     *
c     *  fix pl1 for operations, where dgz is on and dg3 is off, *
c     *    note dg is on if either dgz or dg3 is on              *
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
c        updates made to fix the h2d,h3d files...kac aug 90...
c        updates made to gloopr - call wrth2d before wrtrad (so ctop ok)
c                     to gloopr - send work array to wrth3d
c                     to wrth3d - to write proper layers of heat..
c                                 (in wrtrad)
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
ctune
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
ctune
c
c          wmo aerosols distributions, and b.p.breiglebs surface
c          albedo scheme.            ..............y.h...sep95
c             to gloopr - call insur2 ... input surface albedo and
c                                         aerosols data
c                         call albaer ... compute albedo and aerosols
c                                         distributions
c                         call grlwsw,gaeros ... data initial
#ifdef SW99UPDATE
c        change to opac aerosol data and algorithm (1998)..y.h.oct99
#endif
c
#include <paramodel.h>
#ifndef RSM
#include <comznl.h>
#include <comfgrid.h>
#include <comfphys.h>
#include <comfver.h>
#include <comgrad.h>
#include <radiag.h>
#else
#include <rscomf.rerun.h>
#include <rscomltb.h>
#include <rscommap.h>
#include <rscomgrad.h>
#endif
c
#ifndef RSM
#ifdef MP
#include <commpi.h>
#define LONF22S lonf22p_
#define LONF2S lonf2p_
#define LONLENS lons2
#define LATG2S latg2p_
#define LATGS latgp_
#else
#define LONF22S lonf22_
#define LONF2S lonf2_
#define LONLENS lons2
#define LATG2S latg2_
#define LATGS latg_
#endif
#else
#ifdef MP
#include <comrmpi.h>
#ifdef RSMVECTORIZE
#define LONF22S igrd12p_*jgrd12p_
#define LONF2S igrd12p_*jgrd12p_
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONF22S igrd12p_
#define LONF2S igrd12p_
#define LONLENS lonlen(mype)*2
#endif
#define LATG2S jgrd12p_
#else
#undef RSMVECTORIZE
#define LONF22S igrd12_
#define LONF2S igrd12_
#define LONLENS igrd12_
#define LATG2S jgrd12_
#endif
#endif
#ifdef NFDR
#include <comnfdr.h>
      parameter(twopm=14.)
#endif
      dimension  pllamgr(LONF22S)
      dimension  plphigr(LONF22S)
      dimension     plgr(LONF22S)
      dimension    uugrs(LONF22S,levs_)
      dimension    vvgrs(LONF22S,levs_)
      dimension    tvgrs(LONF22S,levs_)
      dimension    rqgrs(LONF22S,levh_)
      dimension    xxgrs(LONF22S,levs_)
c
c local array
c
      dimension  plamgr(LONF2S)
      dimension  pphigr(LONF2S)
      dimension     pgr(LONF2S)
      dimension    ugrs(LONF2S,levs_)
      dimension    vgrs(LONF2S,levs_)
      dimension    tgrs(LONF2S,levs_)
      dimension    qgrs(LONF2S,levh_)
      dimension    xgrs(LONF2S,levs_)
      dimension    qcis(LONF2S,levs_)
c
#ifdef DG3
#include <comgda.h>
      dimension gda(nwgda,kdgda)
#endif
#ifdef DG
      dimension cldt(LONF2S,levs_),clcv(LONF2S,levs_)
#endif
#ifdef DGZ
      dimension sfcp(LONF2S)
#endif
#ifdef CONDENS
      dimension taucld(LONF2S,levs_),
     +          cldwp(LONF2S,levs_), cldip(LONF2S,levs_)
#endif
c
      dimension workr(LONF2S,levs_)
      dimension w2(LONF2S,levs_)
      dimension tsear(LONF2S),shelgr(LONF2S),
     *          albdoa(LONF2S),
     *          cldary(LONF2S,levs_),
     *          cldsa(LONF2S,4),
     *          cldtot(LONF2S,levs_),
     *          cldcnv(LONF2S,levs_),
     *          mtopa(LONF2S,3),
     *          mbota(LONF2S,3),
     *          swhr(LONF2S,levs_),
     *          hlwr(LONF2S,levs_),
     *          solinr(LONF2S,levs_),
     *          sfnswr(LONF2S),sfdlwr(LONF2S),
#ifdef VIC
     *          sfdswr(LONF2S),sfuswr(LONF2S),
#endif
     *             tsflwr(LONF2S),
     *             tgr(LONF2S)
      dimension    gdfvbr(LONF2S),gdfnbr(LONF2S),
     *             gdfvdr(LONF2S),gdfndr(LONF2S)
#ifdef SW99UPDATE
c
c new aerosol scheme
c
      dimension    rhrh(LONF2S,levs_),
     1             idxc(nxc,LONF2S),
     *             cmix(nxc,LONF2S),
     2             denn(ndn,LONF2S),kprf(LONF2S)
#endif
c
      dtacc = min(dtswav,dtlwav)*3600.
c
#ifndef RSM
      latco=LATGS+1-lat
#endif
      ncloud = ncldg_
      icloud = ntrac_+ 1
c
      do i=1,LONLENS
#ifndef RSM
        pgr(i) = exp(plgr(i))
#else
        pgr(i) = plgr(i)
#endif
#ifndef SCM
        plamgr(i)=pllamgr(i)
        pphigr(i)=plphigr(i)
#endif
      enddo
c....
      do k=1,levs_
       do i=1,LONLENS
        ugrs(i,k)=uugrs(i,k)
        vgrs(i,k)=vvgrs(i,k)
        tgrs(i,k)=tvgrs(i,k)
        xgrs(i,k)=xxgrs(i,k)
       enddo
      enddo
c
      do k=1,levh_
       do i=1,LONLENS
        qgrs(i,k)=rqgrs(i,k)
       enddo
      enddo
c
c    convert virt. temp to thermodynamic temp.
c
      do k=1,levs_
       do i=1,LONLENS
        if(qgrs(i,k).le.0.0) qgrs(i,k)=1.0e-10
        w2(i,k)=1.+0.6*qgrs(i,k)
         qcis(i,k) = 0.
         if(ncloud.gt.0) then
            do ic = icloud,ntotal_
              kc = (ic-1)*levs_ + k
#ifndef NISLQ
              qcis(i,k) = max(qgrs(i,kc),qmin_) + qcis(i,k)
#else
              qcis(i,k) = max(qgrs(i,kc),qmin_) + qcis(i,k)
#endif
            enddo
            w2(i,k)=w2(i,k) - qcis(i,k)
         endif
       enddo
      enddo
      do k=1,levs_
       do i=1,LONLENS
        tgrs(i,k)=tgrs(i,k)/w2(i,k)
       enddo
      enddo
c###
c  after here, qcis is used for calculating cloud water path only.
c  so, qcis is calculated with using only cloud water and ice, 
c  without rain and snow water.
      if(ncloud.gt.0) then
        if (ncloud.le.2) then
          icmx=icloud          ! only cloud water/ice
        else
          icmx=icloud+1        ! only cloud water,ice
        endif
        do k=1,levs_
        do i=1,LONLENS
          qcis(i,k) = 0.
          do ic = icloud,icmx
            kc = (ic-1)*levs_ + k
            qcis(i,k) = max(qgrs(i,kc),0.0e0) + qcis(i,k)
          enddo
        enddo
        enddo
      endif
c###
c
      do i=1,LONLENS
         tsear(i)=tsea(i,lat)
        shelgr(i)=sheleg(i,lat)
           tgr(i)=stc(i,lat,1)
      enddo
c
#ifdef DBG
      print *,'maxmin from grphys at lat=',lat
      call maxmin(tgrs,LONF2S,levs_,1,levs_,'temp')
      call maxmin(qgrs,LONF2S,levh_,1,levh_,'q')
      call maxmin(tsear(1),LONLENS,1,1,1,'tsear')
      call maxmin(shelgr(1),LONLENS,1,1,1,'shelgr')
      call maxmin(tgr(1),LONLENS,1,1,1,'tgr')
#ifdef RSM
      call maxmin(albed (1,lat),LONF2S,1,1,1,'albed')
#endif
      call maxmin(slmsk (1,lat),LONF2S,1,1,1,'slmsk')

#endif
c
      call grrad1(LONLENS,
#ifdef SCM
     &            xgrs(1,1),
#else
     &            xgrs(1,1),plamgr(1),pphigr(1),
#endif
     1            ugrs(1,1),vgrs(1,1),tgrs(1,1),
     2            qgrs(1,1),pgr(1),
#ifndef RSM
     3            albedr(1,lat),slmskr(1,lat),
     4            xlon(1,lat),xlat(1,lat),
#else
     3            albed (1,lat),slmsk (1,lat),
     4            rlon(1,lat),rlat(1,lat),
#endif
     5            tsear(1),shelgr(1),tgr(1),
     6            cvr(1,lat),cvtr(1,lat),cvbr(1,lat),rhcl,
     7            ozon(1,1,lat),albdoa(1),cldary(1,1),
#ifndef SW99UPDATE
     8            cldtot(1,1),cldcnv(1,1),
#else
     8            cldtot(1,1),cldcnv(1,1),rhrh(1,1),
#endif
     9            cldsa(1,1),mtopa(1,1),mbota(1,1),
#ifndef RSM
     a            rbs2,lat,latco,latrue,istrat,
#else
     a            rrs2,lat,latco,latrue,istrat,
#endif
     b            kalb,jo3,slag,rsin1,rcos1,rcos2,
     c            fjd,dlt,jsno,workr,levs_
#ifdef DG
     d            ,cldt(1,1),clcv(1,1)
#endif
#ifdef DGZ
     e            ,sfcp(1)
#endif
#ifdef SW99UPDATE
     b            ,icwp,qcis(1,1)
#endif
     g            )
c
c  use forecasted ozone if possible
c
#ifndef O3C
      if (ntrac_.ge.2.and.thour.gt.0.) then
        do k=1,levs_
          do i=1,LONLENS
#ifdef RSMVECTORIZE
            ozon(i,lat,k) = qgrs(i,k+levs_)
#else
            ozon(i,k,lat) = qgrs(i,k+levs_)
#endif
          enddo
        enddo
      endif
#endif
#ifdef SW99UPDATE
c
c  map global aerosol data to model grid
c
      if (iswsrc(1) .gt. 0) then
        rdg  = 180.0 / pi_
        hdlt = 0.5 * dltgae
c
        i2 = 1
        j2 = 1
        do i=1,LONLENS
#ifndef RSM
          tmp1 = xlon(i,lat) * rdg
#else
          tmp1 = rlon(i,lat) * rdg
#endif
          do i1=1,imxae
            tmp2 = dltgae * (i1 - 1)
            if (tmp2 .gt. 360.0-hdlt) then
              tmp2 = tmp2 - 360.0
            end if
            if (abs(tmp1-tmp2) .le. hdlt) then
              i2 = i1
              go to 40
            end if
          end do
#ifndef RSM
  40      tmp1 = xlat(i,lat) * rdg
#else
  40      tmp1 = rlat(i,lat) * rdg
#endif
          do j1=1,jmxae
            tmp2 = 90.0 - dltgae * (j1 - 1)
            if (abs(tmp1-tmp2) .le. hdlt) then
              j2 = j1
              go to 50
            end if
          end do
c
  50      kprf(i) = kprfg(i2,j2)
          do k=1,nxc
            idxc(k,i) = idxcg(k,i2,j2)
            cmix(k,i) = cmixg(k,i2,j2)
          end do
          do k=1,ndn
            denn(k,i) = denng(k,i2,j2)
          end do
        end do
      end if
c
#endif
      swhr(1:LONF2S,1:levs_)=0.0e0
      hlwr(1:LONF2S,1:levs_)=0.0e0
      solinr(1:LONF2S,1:levs_)=0.0e0
      call grrad2(LONLENS,
     &            tgrs(1,1),qgrs(1,1),qcis(1,1),pgr(1),
#ifndef SW99UPDATE
     1            paerr(1,1,lat),ozon(1,1,lat),albdoa(1),
#else
     1            ozon(1,1,lat),albdoa(1),
#endif
#ifndef RSM
     2            slmskr(1,lat),coszer(1,lat),coszdg(1,lat),
     3            xlat(1,lat),tsear(1),
#else
     2            slmsk (1,lat),coszer(1,lat),coszdg(1,lat),
     3            rlat(1,lat),tsear(1),
#endif
     4            alvbr(1,lat),alnbr(1,lat),alvdr(1,lat),alndr(1,lat),
     5            cldary(1,1),cldtot(1,1),cldcnv(1,1),cldsa(1,1),
     6            mtopa(1,1),mbota(1,1),
#ifdef SW99UPDATE
     6            rhrh(1,1),kprf(1),idxc(1,1),cmix(1,1),denn(1,1),
#endif
#ifndef RSM
     7            lat,latco,sdec,solc,rsin1,rcos1,rcos2,
#else
     7            lat,sdec,solc,rsin1,rcos1,rcos2,
#endif
     8            raddt,dtlw,itimsw,itimlw,kalb,
     8            iswsrc,ibnd,ko3,icwp,icfc,
     9            swhr(1,1),hlwr(1,1),
#ifdef CONDENS
     +            taucld(1,1),
     +            cldwp(1,1), cldip(1,1),
#endif
     a            sfnswr(1),sfdlwr(1),tsflwr(1),
#ifdef VIC
     a            sfdswr(1),sfuswr(1),
#endif
     b            gdfvdr(1),gdfndr(1),
     c            gdfvbr(1),gdfnbr(1),solinr(1,1))
c
#ifdef NFDR
c
c  store 2pm clouds 
c
      dtrhr=dtacc/3600.
      do i=1,LONLENS
#ifdef RSM
        rloctime=solhr+rlon(i,lat)/0.01745329/15.
#else
        rloctime=solhr+xlon(i,lat)/0.01745329/15.
#endif
        if(rloctime.lt. 0.) rloctime=rloctime+24.
        if(rloctime.ge.24.) rloctime=mod(rloctime,24.)
        if(rloctime-dtrhr.le.twopm.and.rloctime+dtrhr.gt.twopm) then
          cloud2pm(i,lat)=cldsa(i,4)
        endif
      enddo
#endif
c
#ifdef DG
         do 351 k=1,levs_
          do 351 i=1,LONLENS
           if (cldt(i,k).lt.0.e0) then
            cldt(i,k) = 0.e0
           else if (cldt(i,k).gt.1.e0) then
            cldt(i,k) = 100.e0
           else
            cldt(i,k) = cldt(i,k) * 100.e0
           end if
  351    continue
         do 352 k=1,levs_
          do 352 i=1,LONLENS
           if (clcv(i,k).lt.0.e0) then
            clcv(i,k) = 0.e0
           else if (clcv(i,k).gt.1.e0) then
            clcv(i,k) = 100.e0
           else
            clcv(i,k) = clcv(i,k) * 100.e0
           end if
  352    continue
#endif
#ifdef DGZ
         call znlacm(LONLENS,LONF2S,cldt(1,1),sfcp(1),dtacc,nmtcld,lat)
         call znlacm(LONLENS,LONF2S,clcv(1,1),sfcp(1),dtacc,nmtccv,lat)
#endif
#ifdef DG3
         call accdia(LONLENS,LONF2S,cldt(1,1),
     1               dtacc,kdtcld,gda(1,1))
         call accdia(LONLENS,LONF2S,clcv(1,1),
     1               dtacc,kdtccv,gda(1,1))
#endif
c
      if(itimsw.eq.1) then
        do k=1,levs_
          do i=1,LONLENS
#ifdef RSMVECTORIZE
            swh(i,lat,k)=swhr(i,k)
            solin(i,lat,k)=solinr(i,k)
#else
            swh(i,k,lat)=swhr(i,k)
            solin(i,k,lat)=solinr(i,k)
#endif
          enddo
        enddo
        do i=1,LONLENS
          sfcnsw(i,lat)=sfnswr(i)
#ifdef VIC
          sfcdsw(i,lat)=sfdswr(i)
          sfcusw(i,lat)=sfuswr(i)
#endif
        enddo
      endif
c
      if(itimlw.eq.1) then
        do k=1,levs_
          do i=1,LONLENS
#ifdef RSMVECTORIZE
            hlw(i,lat,k)=hlwr(i,k)
#else
            hlw(i,k,lat)=hlwr(i,k)
#endif
          enddo
        enddo
        do i=1,LONLENS
          sfcdlw(i,lat)=sfdlwr(i)
           tsflw(i,lat)=tsflwr(i)
        enddo
      endif
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      numlevs=1
      call maxmin(swh(1,1,lat),LONF2S,levs_,1,numlevs,'swh')
      print *,'lwh in grphys'
      call maxmin(hlw(1,1,lat),LONF2S,levs_,1,numlevs,'lwh')
#ifdef MP
      endif
#endif
#endif
c
      do 650 i=1,LONLENS
        cvavg(i,lat) = cvavg(i,lat) + raddt * cvr(i,lat)
  650 continue
c
#ifdef RSMVECTORIZE
      call fillarray(swh,levs_)
      call fillarray(hlw,levs_)
      call fillarray(sfcnsw,1)
#ifdef VIC
      call fillarray(sfcdsw,1)
      call fillarray(sfcusw,1)
#endif
      call fillarray(sfcdlw,1)
      call fillarray(tsflw,1)
      call fillarray(cvavg,1)
#endif
c
      return
      end
