#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#include <define.h>
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
      subroutine progt2(ims2,imx2,kmx,
     &       rhscnpy,rhsmc,ai,bi,ci,smc,slimsk,
     &       canopy,precip,runoff,snowmt,
#ifndef HYDRO
     &       zsoil,soiltyp,sigmaf,delt,lat)
#else
     &       zsoil,soiltyp,sigmaf,delt,lat,hydrow)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(funcdf,funckt,thsat)
#endif
#include <paramodel.h>
      parameter (scanop=2.,rhoh2o=1000.)
      parameter (ctfil1=.5,ctfil2=1.-ctfil1)
      parameter (rffact=.15)
      dimension rhscnpy(imx2),rhsmc(imx2,kmx)
      dimension ai(imx2,kmx),bi(imx2,kmx),ci(imx2,kmx)
      dimension smc(imx2,kmx),canopy(imx2),precip(imx2),soiltyp(imx2)
      dimension sigmaf(imx2),slimsk(imx2),runoff(imx2)
      dimension zsoil(imx2,kmx),snowmt(imx2)
      integer soiltyp
      real inf, infmax, ksat
c local array
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
      logical flag(ILOTS)
      dimension prcp(ILOTS),inf(ILOTS),
     *          infmax(ILOTS)
      dimension tsat(ILOTS),dsat(ILOTS),
     *          ksat(ILOTS)
      dimension smsoil(ILOTS,lsoil_),
     *          cnpy(ILOTS),
     *          dew(ILOTS)
#ifdef HYDRO
      dimension hydrow(ILOTS,lsoil_)
#endif
      im = ims2
      km = kmx
      latd = 44
      lond = 353
      delt2 = delt * 2.
c
c  precipitation rate is needed in unit of kg m-2 s-1
c
      do i = 1, im
        prcp(i) = rhoh2o * (precip(i)/delt+snowmt(i))
        runoff(i) = 0.
        cnpy(i) = canopy(i)
      enddo
c
c  update canopy water content
c
      do i = 1, im
        if(slimsk(i).eq.1.) then
          rhscnpy(i) = rhscnpy(i) + sigmaf(i) * prcp(i)
          canopy(i) = canopy(i) + delt * rhscnpy(i)
          canopy(i) = max(canopy(i),0.)
          prcp(i) = prcp(i) * (1. - sigmaf(i))
          if(canopy(i).gt.scanop) then
            drip = canopy(i) - scanop
            canopy(i) = scanop
            prcp(i) = prcp(i) + drip / delt
          endif
c
c  calculate infiltration rate
c
          inf(i) = prcp(i)
          tsat(i) = thsat(soiltyp(i))
#ifdef HYDRO
          hydrow(i,1)=hydrow(i,1)-prcp(i)
#else
          infmax(i) = (-zsoil(i,1)) *
     &                ((tsat(i) - smc(i,1)) / delt - rhsmc(i,1))
     &                * rhoh2o
          infmax(i) = max(rffact*infmax(i),0.)
          if(inf(i).gt.infmax(i)) then
            runoff(i) = inf(i) - infmax(i)
            inf(i) = infmax(i)
          endif
          inf(i) = inf(i) / rhoh2o
          rhsmc(i,1) = rhsmc(i,1) - inf(i) / zsoil(i,1)
#endif
        endif
      enddo
c
c  we currently ignore the effect of rain on sea ice
c
      do i = 1, im
        flag(i) = slimsk(i).eq.1.
      enddo
c
c  solve the tri-diagonal matrix
c
#ifndef HYDRO
      do k = 1, km
        do i = 1, im
          if(flag(i))  then
            rhsmc(i,k) = rhsmc(i,k) * delt
            ai(i,k) = ai(i,k) * delt
            bi(i,k) = 1. + bi(i,k) * delt
            ci(i,k) = ci(i,k) * delt
          endif
        enddo
      enddo
#else
      do i = 1, im
        if(flag(i))  then
          rhsmc(i,1) = rhsmc(i,1) * delt
          rhsmc(i,2) = rhsmc(i,2) * delt
          ai(i,1) = 1. + ai(i,1) * delt / 
     1                   (-zsoil(i,1)*rhoh2o)
          bi(i,1) = bi(i,1) * delt /
     1                   (-zsoil(i,1)*rhoh2o)
          ci(i,1) = smc(i,1) + (ci(i,1) - hydrow(i,1))*delt / 
     1                   (-zsoil(i,1)*rhoh2o)
          ai(i,2) = ai(i,2) * delt /
     1                 (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
          bi(i,2) = 1. + bi(i,2) * delt / 
     1                 (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
          ci(i,2) = smc(i,2) + (ci(i,2) - hydrow(i,2))*delt / 
     1                 (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
        endif
      enddo
#endif
#ifndef HYDRO
c
c  implicit solver
c
c  forward elimination
c
      do i = 1, im
        if(flag(i)) then
          ci(i,1) = -ci(i,1) / bi(i,1)
          rhsmc(i,1) = rhsmc(i,1) / bi(i,1)
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if(flag(i)) then
            cc = 1. / (bi(i,k) + ai(i,k) * ci(i,k-1))
            ci(i,k) = -ci(i,k) * cc
            rhsmc(i,k) = (rhsmc(i,k) - ai(i,k) * rhsmc(i,k-1)) * cc
          endif
        enddo
      enddo
c  backward substituttion
      do i = 1, im
        if(flag(i)) then
          ci(i,km) = rhsmc(i,km)
        endif
      enddo
      do k = km-1, 1
        do i = 1, im
          if(flag(i)) then
            ci(i,k) = ci(i,k) * ci(i,k+1) + rhsmc(i,k)
          endif
        enddo
      enddo
c
c  update soil moisture
c
      do k = 1, km
        do i = 1, im
          if(flag(i)) then
            smsoil(i,k) = smc(i,k) + ci(i,k)
            smsoil(i,k) = max(smsoil(i,k),0.)
            tdif = max(smsoil(i,k) - tsat(i),0.)
            runoff(i) = runoff(i) -
     &                rhoh2o * tdif * zsoil(i,k) / delt
            smsoil(i,k) = smsoil(i,k) - tdif
          endif
        enddo
      enddo
#else
c
c explicit solver
c
      do i=1,im
        if(flag(i)) then
          if (km.ne.2) then
            print *,'This code is only for two soil levels'
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
c
c acr
c solving the following system of two equations and two unknowns 
c with x1=top layer soil moisture content and x2=bottom layer 
c soil moisture content:
c ai(i,1)*x1 + bi(i,1)*x2 = ci(i,1)
c ai(i,2)*x1 + bi(i,2)*x2 = ci(i,1)
c
          smsoil(i,1) = ci(i,1)/ai(i,1)-bi(i,1)/ai(i,1) *
     1                 (ci(i,2) - ai(i,2)*ci(i,1)/ai(i,1)) / 
     2		       (bi(i,2) - bi(i,1)*ai(i,2)/ai(i,1))
          smsoil(i,2) = (ci(i,2)-ai(i,2)*ci(i,1)/ai(i,1)) / 
     1                 (bi(i,2)-bi(i,1)*ai(i,2)/ai(i,1))
        endif
      enddo
c
c  update soil moisture
c
      do k = 1, km
        do i = 1, im
          if(flag(i)) then
            if(k.gt.1) then
              dz=zsoil(i,k)-zsoil(i,k-1)
            else
              dz=zsoil(i,k)
            endif
c
c  handling negative soil moisture
c
            if(smsoil(i,k).lt.0.) then
              runoff(i)=runoff(i)-smsoil(i,k)*rhoh2o*dz/delt
              rnof=smsoil(i,k)
              smsoil(i,k)=0.
            else
              rnof=0.
            endif
c
c  handling over-saturation of soil moisture
c
            tdif = max(smsoil(i,k) - tsat(i),0.)
            runoff(i) = runoff(i)-rhoh2o*tdif*dz/delt
            smsoil(i,k) = smsoil(i,k) - tdif
          endif
        enddo
      enddo
#endif
      do k = 1, km
        do i = 1, im
          if(flag(i)) then
            smc(i,k) = smsoil(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(flag(i)) then
          canopy(i) = canopy(i)
        endif
      enddo
      return
      end
