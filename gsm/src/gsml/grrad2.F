#define CAL cal
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
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LTSTP ltstp
#define MLVARK mlvark
#define NPES npes
#define RLPNT rlpnt
#define RLTSTP rltstp
#define RMLVARK rmlvark
#define RSLVARK rslvark
#define SLVARK slvark
#include <define.h>
      subroutine grrad2(lons2,
#ifndef SW99UPDATE
     &                  tgrs,qgrs,qcis,pgr,paerr,ozonea,albdoa,
#else
#ifndef RSM
     &                  tgrs,qgrs,qcis,pgr,ozonea,albdoa,
#else
     &                  tgrs,qgrs,qcis,pgr,paerr,ozonea,albdoa,
#endif
#endif
     1                  slmskr,coszro,coszdg,rlat,tsear,
     2                  alvbr,alnbr,alvdr,alndr,
     3                  cldary,cldtot,cldcnv,cldsa,mtopa,mbota,
#ifdef SW99UPDATE
     3                  rhrh,kprf,idxc,cmix,denn,
#endif
#ifndef RSM
     4                  lat,latco,sdec,solc,rsin1,rcos1,rcos2,
#else
     4                  lat,sdec,solc,rsin1,rcos1,rcos2,
#endif
     5                  raddt,dtlw,itimsw,itimlw,kalb,
     5                  iswsrc,ibnd,ko3,icwp,icfc,
     6                  swhr,hlwr,
#ifdef CONDENS
     +                  taucld,cldwp,cldip,
#endif
     6                  sfnswr,sfdlwr,tsflwr,
#ifdef VIC
     6                  sfdswr,sfuswr,
#endif
     7                  gdfvdr,gdfndr,gdfvbr,gdfnbr,solinr)
c     7                  gdfvdr,gdfndr,gdfvbr,gdfnbr)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c     ************************************************************
c     *  added accumulation of clds and convective cloud in dg3  *
c     *                                        k.a.c sept 1994   *
c     *  f3d added for clouds..mi new code=f94/source2/diagnew   *
c     *                               b katz + k.a.c oct  1994   *
c     *  changed kenpts to store total cloud and all lyrs        *
c     *     of cloud......                                       *
c     *                                         k.a.c. nov94     *
c     *  fix pl1 for operations, where dgz is on and dg3 is off, *
c     *     ....note dg is on if either dgz or dg3 is on         *
c     *                                         k.a.c. jan94     *
c     ************************************************************
c
c        updates made to fix the h2d,h3d files...kac aug 90...
c        updates made to gloopr - call wrth2d before wrtrad (so ctop ok)
c                     to gloopr - send work array to wrth3d
c                     to wrth3d - to write proper layers of heat..
c                                 (in wrtrad)
c        updates made to add grid point diagnostics ..k.a.c...sep 91
c                     to gloopr -
c        updates made to pass and receive sib data  ..k.a.c...mar 92
c                     to gloopr -
c        updates made to fix sw rad diagnostics     ..k.a.c...jun 92
c                     proper diurnal weighting
c                     to gloopr and coszmn
c        updates made to calculate clear-sky "on-the-fly" kac aug 92
c                     to gloopr,radfs,fst,spa,lwr,swr
c                     ...for cloud forcing....
c        updates made to call cld optical property routine (cldprp),
c                     to give cld emissivity, optical depth, layer
c                     reflectance and transmitance
c                     to gloopr and radfs ...y.h.           ...feb93
c
c        updates made to change definition of h,m,l domains..
c          to cldprp   k.a.c...dec92 + aug93
c
c          wmo aerosols distributions, and b.p.breiglebs surface
c          albedo scheme.            ..............y.h...sep95
c             to grrad2 - call rdlwsw ... to replace the old radfs
c                                     use gfdls lw and chous sw
c                         kalb=0 use old sfc albedo scheme (matthews)
c                             =1 use new sfc albedo scheme (breigleb)
c                         iswsrc flags for selections of sw absorbers
c                                1:aerosol, 2:o2, 3:co2, 4:w.v. 5:o3
c                                 =1: with the absorber, =0: without
#ifdef SW99UPDATE
c        updates made to opac aerosol algorithm (1998)...y.h..sep99
#endif
c
#include <paramodel.h>
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
#define LONLENS lons2
#include <comfver.h>
#include <radiag.h>
#include <comgpd.h>
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#else
#include <rscomver.h>
#include <rscomrad.h>
#include <rscomgpd.h>
#ifdef MP
#include <comrmpi.h>
#ifdef RSMVECTORIZE
#define LONF2S igrd12p_*jgrd12p_
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONF2S igrd12p_
#define LONLENS lonlen(mype)*2
#endif
#define LATG2S jgrd12p_
#else
#define LONF2S  igrd12_
#define LONLENS igrd12_
#define LATG2S  jgrd12_
#endif
#endif
      parameter (cnwatt=-cal_*1.e4/60.,cnprog=1./cnwatt)
      dimension tgrs(LONF2S,levs_),
     *          qgrs(LONF2S,levs_)
      dimension pgr(LONF2S),
     *          qcis(LONF2S,levs_)
      dimension coszer(LONF2S),coszdg(LONF2S)
      dimension rlat(LONF2S),slmskr(LONF2S)
c
      dimension solinr(LONF2S,levs_) 
#ifndef SW99UPDATE
      dimension cldary(LONF2S,levs_)
      dimension cldtot(LONF2S,levs_),
     *          cldcnv(LONF2S,levs_)
#else
#include <aerparm.h>
      dimension cldary(LONF2S,levs_),
     *          rhrh (LONF2S,levs_)
      dimension cldtot(LONF2S,levs_),
     *          cldcnv(LONF2S,levs_)
      dimension idxc(nxc,LONF2S),
     *          cmix(nxc,LONF2S)
      dimension denn(ndn,LONF2S),
     *          kprf(LONF2S)
#endif
c
      dimension flwup(LONF2S),fswup(LONF2S),
     *          fswdn(LONF2S)
      dimension sswup(LONF2S),sswdn(LONF2S),
     1          slwup(LONF2S),slwdn(LONF2S)
#ifdef CLR
      dimension flwup0(LONF2S),fswup0(LONF2S)
      dimension sswup0(LONF2S),sswdn0(LONF2S),
     *          slwdn0(LONF2S)
#endif
c
#ifndef SW99UPDATE
      dimension    paerr(LONF2S,5),
#else
#ifdef RSM
      dimension    paerr(LONF2S,5)
#endif
      dimension
#endif
     *    ozonea(LONF2S,levs_),
     *    albdoa(LONF2S),
     *    cldsa(LONF2S,4),
     *    mtopa(LONF2S,3),
     *    mbota(LONF2S,3),
     *    coszro(LONF2S),
     *    tsear(LONF2S),
     *    swhr(LONF2S,levs_),
     *    hlwr(LONF2S,levs_),
     *    sfnswr(LONF2S),sfdlwr(LONF2S),
     *    tsflwr(LONF2S)
#ifdef VIC
     *   ,sfdswr(LONF2S),sfuswr(LONF2S)
#endif
      dimension alvbr(LONF2S),alnbr(LONF2S),
     1          alvdr(LONF2S),alndr(LONF2S),
     2          gdfvbr(LONF2S),gdfnbr(LONF2S),
     3          gdfvdr(LONF2S),gdfndr(LONF2S)
#ifdef CONDENS
      dimension    taucld(LONF2S,levs_),
     +             cldwp(LONF2S,levs_),cldip(LONF2S,levs_)
#endif
chybrid
      dimension sihyb(LONF2S,levp1_),cihyb(LONF2S,levp1_)
      dimension slhyb(LONF2S,levs_),clhyb(LONF2S,levs_)
      dimension delhyb(LONF2S,levs_),slkhyb(LONF2S,levs_)
c
c   downward sw fluxes from sw sib rad
c   saved for h2d file
c
#ifdef DG3
#ifdef RSMVECTORIZE
      common/sibsw/ dfvbr(LONF2S,1),
     *              dfnbr(LONF2S,1),
     *              dfvdr(LONF2S,1),
     *              dfndr(LONF2S,1)
#else
      common/sibsw/ dfvbr(LONF2S,LATG2S),
     *              dfnbr(LONF2S,LATG2S),
     *              dfvdr(LONF2S,LATG2S),
     *              dfndr(LONF2S,LATG2S)
#endif
#endif
c
      dimension  iswsrc(nsrc)
#ifdef DBG
c
      if(lat.eq.10) then
      print *, 'before rdlwsw '
      print *,si,sl,lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(10),
     2            qgrs(10,10),tgrs(10,10),ozonea(10,10),tsear(10),
     3            coszro(10),cldary(10,10),albdoa(10),
     4            alvbr(10),alnbr(10),alvdr(10),alndr(10),paerr(10,3)
      print *,'itime ',itimsw,itimlw,ko3,kalb,iswsrc,
     6            hlwr(10,10),slwup(10),slwdn(10),flwup(10),
     7            swhr(10,10),sswup(10),sswdn(10),fswup(10),fswdn(10)
      endif
#endif
c
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
      call rdlwsw(LONLENS,LONF2S,pgr(1),sihyb,slhyb,
     1            lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(1),
     2            qgrs(1,1),tgrs(1,1),ozonea(1,1),tsear(1),
     3            coszro(1),cldary(1,1),cldtot(1,1),cldcnv(1,1),
     4            albdoa(1),
#ifndef SW99UPDATE
     4            alvbr(1),alnbr(1),alvdr(1),alndr(1),paerr(1,1),
#else
     4            alvbr(1),alnbr(1),alvdr(1),alndr(1),
     4            kprf(1),idxc(1,1),cmix(1,1),denn(1,1),rhrh(1,1),
#endif
     5            itimsw,itimlw,ko3,kalb,iswsrc,
     5            ibnd,icwp,qcis,icfc,
     6            hlwr(1,1),slwup,slwdn,flwup,
     7            swhr(1,1),sswup,sswdn,fswup,fswdn,
#ifdef CLR
     8            flwup0,fswup0,sswdn0,sswup0,slwdn0,
#endif
#ifdef CONDENS
     +            taucld,cldwp,cldip,
#endif
     9            gdfvbr(1),gdfnbr(1),gdfvdr(1),gdfndr(1),solinr)
c     9            gdfvbr(1),gdfnbr(1),gdfvdr(1),gdfndr(1))
#ifdef DBG
c
      if(lat.eq.10) then
      print *, 'after rdlwsw '
      print *,si,sl,lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(10),
     2            qgrs(10,10),tgrs(10,10),ozonea(10,10),tsear(10),
     3            coszro(10),cldary(10,10),albdoa(10),
     4            alvbr(10),alnbr(10),alvdr(10),alndr(10),paerr(10,3)
      print *,'itime ',itimsw,itimlw,ko3,kalb,iswsrc,
     6            hlwr(10,10),slwup(10),slwdn(10),flwup(10),
     7            swhr(10,10),sswup(10),sswdn(10),fswup(10),fswdn(10)
      endif
#endif
c
c  cnprog is conversion from w/m**2 to progtm units
c
      if(itimsw.eq.1) then
        do 400 i=1,LONLENS
          sfnswr(i) = (sswdn(i)-sswup(i))*cnprog
#ifdef VIC
          sfdswr(i) = sswdn(i)*cnprog
          sfuswr(i) = sswup(i)*cnprog
#endif
  400   continue
      endif
      if(itimlw.eq.1) then
        do 410 i=1,LONLENS
          sfdlwr(i) = slwdn(i)*cnprog
          tsflwr(i) = tgrs(i,1)
  410   continue
      endif
c
c   save 4 components of downward sw flux
c
#ifdef DG
      if(itimsw.eq.1) then
#endif
c
c  accumulate for h2d file..
c
#ifdef DG
        do 1889 i=1,LONLENS
          dfvbr(i,lat) = dfvbr(i,lat) + raddt * gdfvbr(i)
          dfnbr(i,lat) = dfnbr(i,lat) + raddt * gdfnbr(i)
          dfvdr(i,lat) = dfvdr(i,lat) + raddt * gdfvdr(i)
          dfndr(i,lat) = dfndr(i,lat) + raddt * gdfndr(i)
#endif
 1889   continue
#ifdef DG
      endif
#endif
c
c  grid point monitor-data on radiation grid
c
#ifdef KEN
      if(isave.ne.0.and.npoint.gt.0) then
        do 336 igpt=1,npoint
        if (lat.eq.jgrd(igpt)) then
         do 335 id=1,LONLENS
          if (id.eq.igrd(igpt)) then
           svdata( 25,igpt,itnum)= cldsa(id,3)
           svdata( 26,igpt,itnum)= cldsa(id,2)
           svdata( 27,igpt,itnum)= cldsa(id,1)
           svdata( 38,igpt,itnum)= cldsa(id,4)
           if(isshrt.lt.1.and.itimsw.eq.1) then
            svdata( 41,igpt,itnum)= id
            svdata( 42,igpt,itnum)= lat
            svdata( 43,igpt,itnum)= slmskr(id)
            svdata( 44,igpt,itnum)= tsear (id)
            svdata( 45,igpt,itnum)= sswdn(id)
            svdata( 46,igpt,itnum)= sswup(id)
            svdata( 48,igpt,itnum)= fswdn(id)
            svdata( 49,igpt,itnum)= fswup(id)
            svdata( 51,igpt,itnum)= mtopa(id,3)
            svdata( 52,igpt,itnum)= mtopa(id,2)
            svdata( 53,igpt,itnum)= mtopa(id,1)
            svdata( 54,igpt,itnum)= mbota(id,3)
            svdata( 55,igpt,itnum)= mbota(id,2)
            svdata( 56,igpt,itnum)= mbota(id,1)
            svdata( 57,igpt,itnum)= coszro(id)
            svdata( 58,igpt,itnum)= asin(sdec)*180.e0/3.14159265e0
#ifdef CLR
            svdata( 70,igpt,itnum)= flwup0(id)
            svdata( 71,igpt,itnum)= slwdn0(id)
            if (coszro(id).gt.0.) then
              svdata( 72,igpt,itnum)= fswup0(id)*coszdg(id)/coszro(id)
              svdata( 73,igpt,itnum)= sswdn0(id)*coszdg(id)/coszro(id)
              svdata( 74,igpt,itnum)= sswup0(id)*coszdg(id)/coszro(id)
            else
              svdata( 72,igpt,itnum)=0.
              svdata( 73,igpt,itnum)=0.
              svdata( 74,igpt,itnum)=0.
            endif
#endif
           endif
           if(isshrt.lt.1.and.itimlw.eq.1) then
            svdata( 47,igpt,itnum)= slwdn(id)
            svdata( 50,igpt,itnum)= flwup(id)
           endif
           if (ilshrt.lt.1.and.itimsw.eq.1) then
            do 345 kc=1,levs_
#ifndef SW99UPDATE
             cvcl = float(int(cldary(id,kc))/10)*1.e-3
             if (cvcl.gt.0.0e0) then
              svdata(kc+slvark_+7*levs_,igpt,itnum) = cvcl
             else
              svdata(kc+slvark_+7*levs_,igpt,itnum) =
     1              amod(cldary(id,kc),2.e0)
           end if
#else
            svdata(kc+slvark_+7*levs_,igpt,itnum) =
     1        max(cldcnv(id,kc), cldtot(id,kc))
#endif
#endif
345         continue
#ifdef KEN
           endif
          endif
#endif
335      continue
#ifdef KEN
        endif
#endif
336     continue
#ifdef KEN
      endif
#endif
c
c save gridded radiative fluxes and cld data
c
      if(itimsw.eq.1) then
        do 420 i=1,LONLENS
          if(kalb.eq.0) then
             fluxr(i,lat,17) = fluxr(i,lat,17) + raddt * albdoa(i)
          else
             fluxr(i,lat,17) = fluxr(i,lat,17) + raddt *
     1        (0.466e0 * alvdr(i) + 0.529e0 * alndr(i))
          endif
          fluxr(i,lat,26) = fluxr(i,lat,26) + raddt * cldsa(i,4)
  420   continue
      endif
      if(itimlw.eq.1) then
         do i = 1, LONLENS
           fluxr(i,lat,1 ) = fluxr(i,lat,1 ) + dtlw * flwup(i)
           fluxr(i,lat,19) = fluxr(i,lat,19) + dtlw * slwdn(i)
           fluxr(i,lat,20) = fluxr(i,lat,20) + dtlw * slwup(i)
#ifdef CLR
           fluxr(i,lat,21) = fluxr(i,lat,21) + dtlw * flwup0(i)
           fluxr(i,lat,25) = fluxr(i,lat,25) + dtlw * slwdn0(i)
#endif
         enddo
      endif
c
c   proper diurnal sw wgt..coszro=mean cosz over daylight, while
c                           coszdg= mean cosz over entire interval
      do 3420 i=1,LONLENS
       if (coszro(i).gt.0.) then
        fluxr(i,lat,2 ) = fluxr(i,lat,2 ) + raddt * fswup(i)
     1                                  * coszdg(i)/coszro(i)
        fluxr(i,lat,3 ) = fluxr(i,lat,3 ) + raddt * sswup(i)
     1                                  * coszdg(i)/coszro(i)
        fluxr(i,lat,4 ) = fluxr(i,lat,4 ) + raddt * sswdn(i)
     1                                  * coszdg(i)/coszro(i)
        fluxr(i,lat,18) = fluxr(i,lat,18) + raddt * fswdn(i)
     1                                  * coszdg(i)/coszro(i)
#ifdef CLR
        fluxr(i,lat,22) = fluxr(i,lat,22) + raddt * fswup0(i)
     1                                  * coszdg(i)/coszro(i)
        fluxr(i,lat,23) = fluxr(i,lat,23) + raddt * sswdn0(i)
     1                                  * coszdg(i)/coszro(i)
        fluxr(i,lat,24) = fluxr(i,lat,24) + raddt * sswup0(i)
     1                                  * coszdg(i)/coszro(i)
#endif
       end if
 3420 continue
c
c   save cld frac,toplyr,botlyr and top temp
c   note that order of high, middle and low clouds is
c   reversed for proper output to sflux and h2d files.
c
      do 440 k=1,3
        do 440 i=1,LONLENS
        fluxr(i,lat,8-k) = fluxr(i,lat,8-k) + raddt * cldsa(i,k)
c
c  save interface pressure (cb) of top/bot,
c
        itop = mtopa(i,k)
        ibtc = mbota(i,k)
        fluxr(i,lat,11-k) = fluxr(i,lat,11-k) + raddt *
     1                       sihyb(i,itop+1) * pgr(i)
     2                       * cldsa(i,k)
        fluxr(i,lat,14-k) = fluxr(i,lat,14-k) + raddt *
     1                       sihyb(i,ibtc) * pgr(i)
     2                       * cldsa(i,k)
        fluxr(i,lat,17-k) = fluxr(i,lat,17-k) + raddt *
     1                       tgrs(i,itop) * cldsa(i,k)
  440 continue
c
      return
      end
