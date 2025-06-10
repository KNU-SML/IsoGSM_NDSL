#define CAL cal
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#include <define.h>
      subroutine sfcdrv(ims2,imx2,kmx,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
#ifdef NOAHYDRO
     &               drain,evap,hflx,ep,wind,ec,
#else
     &               drain,evap,hflx,ep,wind,
#endif
#ifdef NOALSM1
     +               snwdph,slc,snoalb,slope,shdmin,shdmax,
#endif
#ifdef VICLSM1
     +               nsl,month,dswflx,binf,ds,dsm,ws,cef,
     +               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,
     +               wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,
     +               apnld,btnld,gmnld,flaild,vrtld,lstsnld,
     +               silz,snwz,csno,rsno,tsf,tpk,sfw,pkw,
#endif
#ifdef VICLSM2
     +               nsl,msub,month,dswflx,binf,ds,dsm,ws,cef,
     +               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,
     +               wpwld,smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,
     +               btnld,gmnld,
     +               nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,
     +               rsnld,tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,
     +               smcld, sicld,stcld,silz,snwz,
#endif
     +               snowfl, runoff, precip, srflag)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
c
#include <paramodel.h>
       integer soiltyp, vegtype, slope
       integer i, j, k
c passing array
      dimension ps(imx2),t1(imx2),q1(imx2)
      dimension sheleg(imx2),snowmt(imx2),snowev(imx2)
      dimension cm(imx2),ch(imx2)
      dimension tskin(imx2),qsurf(imx2),dm(imx2)
      dimension dlwflx(imx2),slrad(imx2)
      dimension smc(imx2,kmx),tg3(imx2),canopy(imx2)
#ifndef VIC
      dimension stc(imx2,kmx)
#else
      dimension stc(imx2,nsl),dswflx(imx2)
#endif
c
      dimension z0rl(imx2),soiltyp(imx2),gflux(imx2),slimsk(imx2)
      dimension drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
      dimension evap(imx2),hflx(imx2),ep(imx2),wind(imx2)
      dimension vegtype(imx2),runoff(imx2)
      dimension precip(imx2), srflag(imx2),snowfl(imx2)
#ifdef NOALSM1
      dimension snoalb(imx2), slope(imx2), shdmin(imx2), shdmax(imx2)
      dimension slc(imx2,kmx), snwdph(imx2)
#endif
#ifdef VICLSM1
      dimension binf(imx2),ds(imx2),dsm(imx2),ws(imx2),cef(imx2)
      dimension expld(imx2,kmx),kstld(imx2,kmx),
     1          dphld(imx2,kmx),bubld(imx2,kmx),
     2          qrtld(imx2,kmx),bkdld(imx2,kmx),
     3          sldld(imx2,kmx),wcrld(imx2,kmx),
     4          wpwld(imx2,kmx),smrld(imx2,kmx),
     5          smxld(imx2,kmx),sicld(imx2,kmx),
     6          dpnld(imx2,nsl),sxnld(imx2,nsl),
     7          epnld(imx2,nsl),bbnld(imx2,nsl),
     8          apnld(imx2,nsl),btnld(imx2,nsl),
     9          gmnld(imx2,nsl)
      real      kstld(imx2,kmx)
      dimension flaild(imx2), vrtld(imx2,kmx),lstsnld(imx2)
      integer lstsnld
      dimension silz(imx2),snwz(imx2),csno(imx2),rsno(imx2)
      dimension tsf(imx2),tpk(imx2),sfw(imx2),pkw(imx2)
#endif
#ifdef VICLSM2
      dimension binf(imx2),ds(imx2),dsm(imx2),ws(imx2),cef(imx2)
      dimension expld(imx2,kmx),
     1          dphld(imx2,kmx),bubld(imx2,kmx),
     2          qrtld(imx2,kmx),bkdld(imx2,kmx),
     3          sldld(imx2,kmx),wcrld(imx2,kmx),
     4          wpwld(imx2,kmx),smrld(imx2,kmx),
     5          smxld(imx2,kmx),
     6          dpnld(imx2,nsl),sxnld(imx2,nsl),
     7          epnld(imx2,nsl),bbnld(imx2,nsl),
     8          apnld(imx2,nsl),btnld(imx2,nsl),
     9          gmnld(imx2,nsl)
      real      kstld(imx2,kmx)
c
      dimension nvegld(imx2), flaild(imx2,msub)
      dimension vfrld(imx2,msub),vtypld(imx2,msub)
      dimension cnpld(imx2,msub),snold(imx2,msub),csnld(imx2,msub)
      dimension rsnld(imx2,msub),tsfld(imx2,msub),tpkld(imx2,msub)
      dimension sfwld(imx2,msub),pkwld(imx2,msub),lstsnld(imx2,msub)
      dimension vrtld(imx2,kmx,msub)
      dimension smcld(imx2,kmx,msub)
      dimension sicld(imx2,kmx,msub)
      dimension stcld(imx2,nsl,msub)
      dimension silz(imx2),snwz(imx2)
c
      integer nvegld,vtypld,lstsnld
#endif
c
c local array
c
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
c
#ifdef NOALSM1
c
c  declare zsoil(4)
c
      dimension zsoil_noah(4)
      data zsoil_noah/-0.1, -0.4, -1.0, -2.0/
#endif
c
#ifdef DBG
      character*4 sfcty
#ifdef OSULSM1
         sfcty = 'osu1'
#endif
#ifdef OSULSM2
         sfcty = 'osu2'
#endif
#ifdef NOALSM1
         sfcty = 'noa1'
#endif
#ifdef VICLSM1
         sfcty = 'vic1'
#endif
#ifdef VICLSM2
         sfcty = 'vic2'
#endif
#endif
c
      im = ims2
      km = kmx
c
c initialize
c
#ifdef DBG
      print *,'--- enter sfcdrv ---'
#endif      

      do i = 1, im
         runoff(i) = 0.
         drain(i) = 0.
         dm(i) = 0.
         evap(i) = 0.
         hflx(i) = 0.
         ep(i) = 0.
         snowfl(i) = 0.
         snowmt(i) = 0.
         snowev(i) = 0.
         gflux(i) = 0.
      enddo
c
c set up soil layer configuration
c
      do i = 1, im
        if(slimsk(i).eq.0.) then
          zsoil(i,1) = 0.
        elseif(slimsk(i).eq.1.) then
          zsoil(i,1) = -.10
        else
          zsoil(i,1) = -3. / km
        endif
      enddo
c
#ifdef DBG
      print *,'in sfcdrv before compute zsoil'
#endif

      do k = 2, km
        do i = 1, im
          if(slimsk(i).eq.0.) then
            zsoil(i,k) = 0.
          elseif(slimsk(i).eq.1.) then
#ifdef OSU
c
c  use default soil layer configuration
c
            zsoil(i,k) = zsoil(i,k-1)
     &                   + (-2. - zsoil(i,1)) / (km - 1)
#endif
#ifdef NOALSM1
c
c  override soil layer configuration using noah_zsoil
c
            zsoil(i,k) = zsoil_noah(k)
#endif
#ifdef VIC
            zsoil(i,k) = -dphld(i,k)
#endif
          else
            zsoil(i,k) = - 3. * float(k) / float(km)
          endif
        enddo
      enddo
c
c  surface energy/water balance over land and seaice
c
#ifdef DBG
      print 88, sfcty
 88   format(' in sfcdrv surface type: ',a4)
#endif
c
#ifdef OSU
       call sfc_osu(ims2,imx2,kmx,
     &               ps, t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,          
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
     &               drain,evap,hflx,ep,wind,
     +               snowfl, runoff, precip, srflag)
#endif
c
#ifdef NOALSM1
c
#ifdef DBG
      print*,' in sfcdrv before sfc_noah'
#endif
c
c  call Noah LSM over land (use sfc_noah_ver2.F)
c
       call sfc_noah(ims2,imx2,kmx,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
     &               drain,evap,hflx, ep, wind,
     +               snwdph, slc, snoalb, slope, shdmin, shdmax,
#ifndef NOAHYDRO
     +               snowfl, runoff, precip, srflag)
#else
     +               snowfl, runoff, precip, srflag, ec)
#endif
c
#endif
c
#ifdef VIC
       do i=1,imx2
         dswflx(i)=-dswflx(i)*cal_*1.e04/60.
       enddo
#endif
#ifdef VICLSM1
c
c  call UW VIC over land (use fortran-version VIC)
c
#ifdef DBG
      print*,' in sfcdrv before sfc_vic1'
#endif
c
       call sfc_vic1(ims2,imx2,kmx,lat,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,sigmaf,vegtype,canopy,
     &               dlwflx,dswflx,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,sl1,slk1,slimsk,
     &               drain,evap,hflx, ep, wind,
     +               nsl,month,binf,ds,dsm,ws,cef,
     +               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,
     +               wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,
     +               apnld,btnld,gmnld,flaild,vrtld,lstsnld,
     +               silz,snwz,csno,rsno,tsf,tpk,sfw,pkw,
     +               snowfl, runoff, precip)
c
#endif
c
#ifdef VICLSM2
c
c  call UW VIC over land (use fortran-version VIC)
c
#ifdef DBG
      print*,' in sfcdrv before sfc_vic2'
#endif
c
       call sfc_vic2(ims2,imx2,kmx,lat,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,sigmaf,vegtype,canopy,
     &               dlwflx,dswflx,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,sl1,slk1,slimsk,
     &               drain,evap,hflx, ep, wind,
     +               nsl,msub,month,binf,ds,dsm,ws,cef,
     +               expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,
     +               wpwld,smrld,smxld,dpnld,sxnld,epnld,bbnld,apnld,
     +               btnld,gmnld,
     +               nvegld,flaild,vfrld,vtypld,cnpld,snold,csnld,
     +               rsnld,tsfld,tpkld,sfwld,pkwld,lstsnld,vrtld,
     +               smcld, sicld,stcld,silz,snwz,
     +               snowfl, runoff, precip)
c
#endif
c
c  call OSU LSM over seaice
c
#ifndef OSU
c
#ifdef DBG
      print*,' in sfcdrv before sfc_seaice'
#endif
c
       call sfc_seaice(ims2,imx2,kmx,
     &               ps, t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
     &               drain,evap,hflx,ep,wind,
     &               snowfl, runoff, precip, srflag)
c
#endif
c
c  surface energy/water balance over ocean
c
#ifdef DBG
      print*,' in sfcdrv before sfc_ocean'
#endif
c
       call sfc_ocean(ims2,imx2,kmx,
     &                ps,t1,q1,tskin,qsurf,
     &                dm, cm, ch,
     &                sl1,slk1,slimsk,inistp,lat,
     &                evap,hflx,wind)
c
      return
      end
