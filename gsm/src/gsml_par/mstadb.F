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
#define RD rd
#define RV rv
      subroutine mstadb(imx2,kmx,k1,k2,sl,slk,ps,tenv,qenv,
     &                  klcl,kbot,ktop,tcld,qcld)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpkap,ftdp,ftlcl,fthe,ftma)
#endif
c     subprogram documentation block
c
c subprogram: mstadb       compute moist adiabatic cloud soundings
c
c   author: n phillips    date: nov 1983
c
c abstract: atmospheric columns of temperature and specific humidity
c   are examined by this routine for conditional instability.
c   the test parcel is chosen from the layer between layers k1 and k2
c   that has the warmest potential wet-bulb temperature.
c   excess cloud temperatures and specific humidities are returned
c   where the lifted parcel is found to be buoyant.
c   fast inlinable functions are invoked to compute
c   dewpoint and lifting condensation level temperatures,
c   equivalent potential temperature at the lcl, and
c   temperature and specific humidity of the ascending parcel.
c
c program history log:
c   83-11     phillips
c   91-05-07  iredell             arguments changed, code tidied
c   02-09     yifeng cui added rsm mp    org:sdsc
c
c
c usage:    call mstadb(im,km,k1,k2,sl,slk,ps,tenv,qenv,
c    &                  klcl,kbot,ktop,tcld,qcld)
c
c   input argument list:
c     im       - integer number of atmospheric columns
c     km       - integer number of sigma levels in a column
c     k1       - integer lowest level from which a parcel can originate
c     k2       - integer highest level from which a parcel can originate
c     sl       - real (km) sigma values
c     slk      - real (km) sigma values to the kappa
c     ps       - real (im) surface pressure in kilopascals (cb)
c     tenv     - real (im,km) environment temperatures
c     qenv     - real (im,km) environment specific humidities
c
c   output argument list:
c     klcl     - integer (im) level just above lcl (km+1 if no lcl)
c     kbot     - integer (im) level just above cloud bottom
c     ktop     - integer (im) level just below cloud top
c              - note that kbot(i) gt ktop(i) if no cloud.
c     tcld     - real (im,km) of excess cloud temperatures.
c                (parcel t minus environ t, or 0. where no cloud)
c     qcld     - real (im,km) of excess cloud specific humidities.
c                (parcel q minus environ q, or 0. where no cloud)
c
c subprograms called:
c     fpkap    - function to compute pressure to the kappa power
c     ftdp     - function to compute dewpoint temperature
c     ftlcl    - function to compute lcl temperature
c     fthe     - function to compute equivalent potential temperature
c     ftma     - function to compute parcel temperature and humidity
c
c remarks: all functions are inlined by fpp.
c          nonstandard automatic arrays are used.
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
#include <paramodel.h>
c
      dimension sl(imx2,kmx),slk(imx2,kmx),ps(imx2),
     $   tenv(imx2,kmx),qenv(imx2,kmx),
     &   klcl(imx2),kbot(imx2),ktop(imx2),tcld(imx2,kmx),qcld(imx2,kmx)
c
c  physical parameters
c
      parameter(rd=rd_,rv=rv_)
      parameter(eps=rd/rv,epsm1=rd/rv-1.,ftv=rv/rd-1.)
c
c  local arrays
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
      dimension psk(ILOTS),slkma(ILOTS),
     *          thema(ILOTS)
c
      parameter(nx=151,ny=121)
      dimension tbtma(nx,ny),tbqma(nx,ny)
      common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma
c
c  determine warmest potential wet-bulb temperature between k1 and k2.
c  compute its lifting condensation level.
c
      do i=1,imx2
        psk(i)=fpkap(ps(i))
        slkma(i)=0.
        thema(i)=0.
      enddo
      do k=k1,k2
        do i=1,imx2
          pv=sl(i,k)*ps(i)*qenv(i,k)/(eps-epsm1*qenv(i,k))
          tdpd=tenv(i,k)-ftdp(pv)
          if(tdpd.gt.0.) then
            tlcl=ftlcl(tenv(i,k),tdpd)
            slklcl=slk(i,k)*tlcl/tenv(i,k)
          else
            tlcl=tenv(i,k)
            slklcl=slk(i,k)
          endif
          thelcl=fthe(tlcl,slklcl*psk(i))
          if(thelcl.gt.thema(i)) then
            slkma(i)=slklcl
            thema(i)=thelcl
          endif
        enddo
      enddo
c-----------------------------------------------------------------------
c  set cloud temperatures and humidities wherever the parcel lifted up
c  the moist adiabat is buoyant with respect to the environment.
      do i=1,imx2
        klcl(i)=kmx+1
        kbot(i)=kmx+1
        ktop(i)=0
      enddo
      do k=1,kmx
        do i=1,imx2
          tcld(i,k)=0.
          qcld(i,k)=0.
        enddo
      enddo
      do k=k1,kmx
        do i=1,imx2
          if(slk(i,k).le.slkma(i)) then
            klcl(i)=min(klcl(i),k)
c
c insert ftma
c           tma=ftma(thema(i),slk(k)*psk(i),qma)
c
            xj=min(max(c1xma+c2xma*thema(i),1.),float(nx))
            yj=min(max(c1yma+c2yma*slk(i,k)*psk(i),1.),float(ny))
            jx=min(xj,nx-1.)
            jy=min(yj,ny-1.)
            ftx1=tbtma(jx,jy)+(xj-jx)*
     1                  (tbtma(jx+1,jy)-tbtma(jx,jy))
            ftx2=tbtma(jx,jy+1)+(xj-jx)*
     2                  (tbtma(jx+1,jy+1)-tbtma(jx,jy+1))
            ftma1=ftx1+(yj-jy)*(ftx2-ftx1)
            qx1=tbqma(jx,jy)+(xj-jx)*(tbqma(jx+1,jy)-tbqma(jx,jy))
            qx2=tbqma(jx,jy+1)+(xj-jx)*(tbqma(jx+1,jy+1)-tbqma(jx,jy+1))
            qma=qx1+(yj-jy)*(qx2-qx1)
            tma=ftma1
c
            tvcld=tma*(1.+ftv*qma)
            tvenv=tenv(i,k)*(1.+ftv*qenv(i,k))
            if(tvcld.gt.tvenv) then
              kbot(i)=min(kbot(i),k)
              ktop(i)=max(ktop(i),k)
              tcld(i,k)=tma-tenv(i,k)
              qcld(i,k)=qma-qenv(i,k)
            endif
          endif
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end
