#define G g
#define IGRD1 igrd1
#define IGRD igrd
#define KO ko
#define LEVS levs
#define RD rd
#define RV rv
      subroutine sig2p(im,ix,km,si,sl,
     &                 ps,us,vs,os,zs,zi,ts,rs,qs,as,cls,prs,o3s,
     &                 ko,po,lpcl,lppr,lpo3,up,vp,op,zp,tp,rp,
     &                 qp,ap,clp,prp,o3p)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c$$$  subprogram documentation block
c
c  subprogram:    sig2p       sigma to pressure interpolation
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: interpolates winds, omega, height, temperature and humidity
c   from the sigma coordinate system to the mandatory pressure levels.
c   assumes that relative humidity, temperature, geopotential heights,
c   wind components and vertical velocity vary linearly in the vertical
c   with the log of pressure.  underground heights are obtained using
c   the shuell method and underground temperatures are obtained using
c   a constant moist adiabatic lapse rate.  heights above the top sigma
c   level are integrated hydrostatically.  otherwise fields are held
c   constant outside the sigma structure and no extrapolation is done.
c
c program history log:
c   92-10-31  sela,newell,gerrity,ballish,deaven,iredell
c
c usage:    call sig2p(im,ix,km,si,sl,
c    &                 ps,us,vs,os,zs,zi,ts,rs,qs,
c    &                 ko,po,up,vp,op,zp,tp,rp)
c
c   input argument list:
c     im       - integer number of points
c     ix       - integer first dimension of upper air data
c     km       - integer number of sigma levels
c     si       - real (km+1) sigma interface values
c     sl       - real (km) sigma values
c     ps       - real (im) surface pressure in kpa
c     us       - real (ix,km) zonal wind in m/s
c     vs       - real (ix,km) merid wind in m/s
c     os       - real (im,km) vertical velocity in pa/s
c     zs       - real (im,km) heights on the full levels in m
c     zi       - real (im,km) heights on the interfaces in m
c     ts       - real (ix,km) temperature in k
c     rs       - real (im,km) relative humidity in percent
c     qs       - real (ix,km) specific humidity in kg/kg
c     ko       - integer number of pressure levels
c     po       - real (ko) mandatory pressures in kpa
c
c   output argument list:
c     up       - real (im,ko) zonal wind in m/s
c     vp       - real (im,ko) merid wind in m/s
c     op       - real (im,ko) vertical velocity in pa/s
c     zp       - real (im,ko) heights in m
c     tp       - real (im,ko) temperature in k
c     rp       - real (im,ko) relative humidity in percent
c
c subprograms called:
c   isrchflt - find first value in an array less than target value
c
c attributes:
c   language: cray fortran
c
c$$$
#include <paramodel.h>
      dimension si(im,km+1),sl(im,km),ps(im)
      dimension us(ix,km),vs(ix,km),os(im,km)
      dimension zs(im,km),zi(im,km),ts(ix,km),rs(im,km),qs(ix,km)
      dimension as(ix,km),o3s(ix,km),cls(ix,km),prs(ix,km)
      dimension po(ko)
c
      dimension up(im,ko),vp(im,ko),op(im,ko)
      dimension zp(im,ko),tp(im,ko),rp(im,ko)
      dimension qp(im,ko),ap(im,ko),o3p(im,ko),clp(im,ko),prp(im,ko)
c
      dimension asi(levs_),asl(levs_),apo(ko_),aps(im)
      parameter(g= g_ ,rd= rd_ ,rv= rv_ )
      parameter(rog=rd/g,fvirt=rv/rd-1.)
      parameter(gammam=-6.5e-3,zshul=75.,tvshul=290.66)
      logical lpo3, lpcl, lppr
c
c  compute log pressures for interpolation
c
      do kk=1,ko
      apo(kk)=log(po(kk))
      enddo
      do ii=1,im
      aps(ii)=log(ps(ii))
      enddo
      apsmin=aps(1)
      apsmax=aps(1)
      do ii=1,im
      apsmin=min(apsmin,aps(ii))
      apsmax=max(apsmax,aps(ii))
      enddo
c 
c  determine sigma layers bracketing pressure layer.
c  within sigma structure, interpolate fields linearly in log pressure
c  between bracketing full sigma layers except heights are interpolated
c  between the nearest full sigma layer and the nearest sigma interface
c
      do i=1,im
        do kk=2,km
          asi(kk)=log(si(i,kk))
        enddo
        do kk=1,km
          asl(kk)=log(sl(i,kk))
        enddo
        kd1=1
        do k=1,ko
          kd1=kd1+isrchflt(km-kd1-1,asl(kd1+1),1,apo(k)-apsmin)-1
          kd2=kd1+isrchflt(km-kd1-1,asl(kd1+1),1,apo(k)-apsmax)-1
          do kd=kd1,kd2
            ku=kd+1
            ask=apo(k)-aps(i)
            if(ask.le.asl(kd).and.ask.gt.asl(ku)) then
              wu=(asl(kd)-ask)/(asl(kd)-asl(kd+1))
              wd=1.-wu
              up(i,k)=wu*us(i,ku)+wd*us(i,kd)
              vp(i,k)=wu*vs(i,ku)+wd*vs(i,kd)
              op(i,k)=wu*os(i,ku)+wd*os(i,kd)
              ki=kd+1
              di=asi(ki)-ask
              kl=nint(ki-0.5+sign(0.5,di))
              wl=di/(asi(ki)-asl(kl))
              wi=1.-wl
              zp(i,k)=wi*zi(i,ki)+wl*zs(i,kl)
              tp(i,k)=wu*ts(i,ku)+wd*ts(i,kd)
              rp(i,k)=wu*rs(i,ku)+wd*rs(i,kd)
              qp(i,k)=wu*qs(i,ku)+wd*qs(i,kd)
              ap(i,k)=wu*as(i,ku)+wd*as(i,kd)
              if(lpcl) clp(i,k)=wu*cls(i,ku)+wd*cls(i,kd)
              if(lppr) prp(i,k)=wu*prs(i,ku)+wd*prs(i,kd)
              if(lpo3) o3p(i,k)=wu*o3s(i,ku)+wd*o3s(i,kd)
            endif
          enddo
        enddo
      enddo
c 
c  interpolate sigma to pressure outside the model domain
c
      do i=1,im
        do kk=2,km
          asi(kk)=log(si(i,kk))
        enddo
        do kk=1,km
          asl(kk)=log(sl(i,kk))
        enddo
        do k=1,ko
          ask=apo(k)-aps(i)
c
c  below ground use shuell method to obtain height, constant lapse rate
c  to obtain temperature, and hold other fields constant
c
          if(ask.gt.0.) then
            up(i,k)=us(i,1)
            vp(i,k)=vs(i,1)
            op(i,k)=os(i,1)
            tvsf=ts(i,1)*(1.+fvirt*qs(i,1))-gammam*(zs(i,1)-zi(i,1))
            if(zi(i,1).gt.zshul) then
              tvsl=tvsf-gammam*zi(i,1)
              if(tvsl.gt.tvshul) then
                if(tvsf.gt.tvshul) then
                  tvsl=tvshul-5.e-3*(tvsf-tvshul)**2
                else
                  tvsl=tvshul
                endif
              endif
              gammas=(tvsf-tvsl)/zi(i,1)
            else
              gammas=0.
            endif
            part=rog*ask
            zp(i,k)=zi(i,1)-tvsf*part/(1.+0.5*gammas*part)
            tp(i,k)=ts(i,1)+gammam*(zp(i,k)-zs(i,1))
            rp(i,k)=rs(i,1)
            qp(i,k)=qs(i,1)
            ap(i,k)=as(i,1)
            if(lpcl) clp(i,k)=cls(i,1)
            if(lppr) prp(i,k)=prs(i,1)
            if(lpo3) o3p(i,k)=o3s(i,1)
c
c  between bottom sigma and ground interpolate height,
c  extrapolate temperature and hold other fields constant
c
          elseif(ask.ge.asl(1)) then
            up(i,k)=us(i,1)
            vp(i,k)=vs(i,1)
            op(i,k)=os(i,1)
            wl=ask/asl(1)
            wi=1.-wl
            zp(i,k)=wi*zi(i,1)+wl*zs(i,1)
            wu=(asl(1)-ask)/(asl(1)-asl(2))
            wd=1.-wu
            tp(i,k)=wu*ts(i,2)+wd*ts(i,1)
            rp(i,k)=rs(i,1)
            qp(i,k)=qs(i,1)
            ap(i,k)=as(i,1)
            if(lpcl) clp(i,k)=cls(i,1)
            if(lppr) prp(i,k)=prs(i,1)
            if(lpo3) o3p(i,k)=o3s(i,1)
c
c  above top sigma integrate height hydrostatically
c  and hold other fields constant
c
          elseif(ask.le.asl(km)) then
            up(i,k)=us(i,km)
            vp(i,k)=vs(i,km)
            op(i,k)=os(i,km)
            tvkm=ts(i,km)*(1.+fvirt*qs(i,km))
            zp(i,k)=zs(i,km)+rog*tvkm*(asl(km)-ask)
            tp(i,k)=ts(i,km)
            rp(i,k)=rs(i,km)
            qp(i,k)=qs(i,km)
            ap(i,k)=as(i,km)
            if(lpcl) clp(i,k)=cls(i,km)
            if(lppr) prp(i,k)=prs(i,km)
            if(lpo3) o3p(i,k)=o3s(i,km)
          endif
        enddo
      enddo
c 
      return
      end
