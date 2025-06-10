#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#include <define.h>

      subroutine sfc_seaice(ims2,imx2,kmx,
     &               ps,t1,q1,sheleg,tskin,qsurf,
     &               smc,stc,dm,soiltyp,sigmaf,vegtype,canopy,
     &               dlwflx,slrad,snowmt,snowev,delt,z0rl,tg3,
     &               gflux,zsoil,cm,ch,
     &               rcl,sl1,slk1,slimsk,inistp,lat,
     &               drain,evap,hflx,ep,wind,
     +               snowfl, runoff,precip,srflag)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0,funcdf,funckt,ktsoil,twlt,thsat)
#endif
c
c ca is the von karman constant
c
#include <paramodel.h>
      integer soiltyp, vegtype
      integer i, j, k
      dimension ps(imx2),                  t1(imx2),q1(imx2)
      dimension sheleg(imx2),snowmt(imx2),snowev(imx2)
      dimension cm(imx2),ch(imx2)
      dimension tskin(imx2),qsurf(imx2),dm(imx2),slrad(imx2)
      dimension smc(imx2,kmx),stc(imx2,kmx),tg3(imx2),canopy(imx2)
      dimension z0rl(imx2),plantr(imx2),soiltyp(imx2),gflux(imx2)
      dimension slimsk(imx2),rhscnpy(imx2),rhsmc(imx2,kmx),rb(imx2)
      dimension aim(imx2,kmx),bim(imx2,kmx),cim(imx2,kmx)
      dimension drain(imx2),zsoil(imx2,kmx),sigmaf(imx2)
      dimension evap(imx2),hflx(imx2),ep(imx2)
      dimension wind(imx2)
      dimension dlwflx(imx2)
      dimension vegtype(imx2)
      dimension runoff(imx2)
      dimension precip(imx2),   srflag(imx2)
      dimension snowfl(imx2)
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
      im = ims2
      km = kmx
c
c  surface energy/water balance over land and seaice
c
c perform snow-rain detection
c
      do 410 j=1, im
      if(slimsk(j) .eq. 2.) then
c
c perform snow-rain detection using precip and srflag
c
       if(srflag(j) .eq. 1.) then
         snowfl(j)=precip(j)
         if(slimsk(j) .ne. 0.) then
           sheleg(j) = sheleg(j) + 1.e3 * precip(j)
         endif
         precip(j) = 0.
      else
         snowfl(j)=0.
      endif
      endif
410   continue
c
c compute sfc energy balance
c
       call seaicetm(ims2, imx2,  kmx,   ps,
     1                     t1,       q1,
     2 sheleg       ,tskin,      qsurf,
     3 smc,   stc,   dm,   soiltyp,sigmaf,vegtype,canopy,
     4 dlwflx,slrad,snowmt,snowev,
     5 delt  ,zorl,       tg3,
     6 gflux,
     7 zsoil,cm,ch,    rhscnpy,rhsmc,aim,bim,cim,
     8 rcl,sl1,slk1,slimsk     ,inistp,lat,
     9 drain,evap,hflx,     ep,              wind)
c
      return
      end
