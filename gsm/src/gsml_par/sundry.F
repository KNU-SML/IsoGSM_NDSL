#undef sun
      subroutine sundry(im,ix,km,kslp,kli,lpcl,lppr,
     &                  clat,si,sl,kt,pt,
     &                  zs,ps,os,u,v,o,t,r,q,qs,tpt,qpt,zm,tm,
     &                  cl,pr,sun)
c
c$$$  subprogram documentation block
c
c subprogram:    sundry      compute sundry fields
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: computes sundry fields.
c   the current nsun=29 sundry fields are:
c     1) surface pressure
c     2) surface pressure tendency
c     3) column precipitable water
c     4) column relative humidity
c     5) tropopause temperature
c     6) tropopause pressure
c     7) tropopause zonal wind
c     8) tropopause meridional wind
c     9) tropopause vertical wind speed shear
c    10) surface lifted index
c    11) best lifted index
c    12) maximum wind level temperature
c    13) maximum wind level pressure
c    14) maximum wind level zonal wind
c    15) maximum wind level meridional wind
c    16) surface orography
c    17) sea level pressure
c
c    18) relative humidity in sigma range (0.44,1.00)
c    19) relative humidity in sigma range (0.72,0.94)
c    20) relative humidity in sigma range (0.44,0.72)
c    21) potential temperature at sigma 0.9950
c    22) temperature at sigma 0.9950
c    23) pressure vertical velocity at sigma 0.9950
c    24) relative humidity at sigma 0.9950
c    25) zonal wind at sigma 0.9950
c    26) meridional wind at sigma 0.9950
c    27) specific humidity at sigma 0.9950
c    28) total cloud water
c    29) total ozone
c
c subprograms called:
c   sig2tp       interpolate sigma to tropopause level
c   sig2mw       interpolate sigma to maxwind level
c   liftix       compute best lifted index
c
c program history log:
c   92-10-31  mccalla,iredell
c
c usage:    call sundry(im,ix,km,kslp,kli,clat,si,sl,kt,pt,
c    &                  zs,ps,os,u,v,o,t,r,q,qs,tpt,qpt,zm,tm,sun)
c
c   input argument list:
c     im       - integer number of points
c     ix       - integer first dimension of upper air data
c     km       - integer number of levels
c     kslp     - integer (2) relevant pressure levels for slp
c     kli      - integer relevant pressure level for lifted index
c     clat     - real (im) cosine of latitude
c     si       - real (km) sigma interfaces
c     sl       - real (km) sigma values
c     kt       - integer number of pressure thickness layers
c     pt       - real pressure thickness in kpa
c     zs       - real (im) surface orography in m
c     ps       - real (im) surface pressure in kpa
c     os       - real (im) surface pressure tendency in pa/s
c     u        - real (ix,km) zonal wind in m/s
c     v        - real (ix,km) meridional wind in m/s
c     o        - real (im,km) vertical velocity in pa/s
c     t        - real (ix,km) temperature in k
c     r        - real (im,km) relative humidity in percent
c     q        - real (ix,km) specific humidity in kg/kg
c     qs       - real (im,km) saturated specific humidity in kg/kg
c     tpt      - real (im,kt) temperature in k
c     qpt      - real (im,kt) specific humidity in kg/kg
c     zm       - real (im,*) height on pressure surface in m
c     tm       - real (im,*) temperature on pressure surface in k
c     cl       - real (ix,*) cloud water in kg/kg
c     pr       - real (ix,*) precipitation water in kg/kg
c
c   output argument list:
c     sun      - real (im,nsun) sundry fields given above
c
c attributes:
c   language: cray fortran
c
c$$$
      parameter(nsun=29)
      parameter(nps  = 1,npst = 2,ntpw = 3,ntrh = 4,ntpt = 5,
     &          ntpp = 6,ntpu = 7,ntpv = 8,ntpsh= 9,nsli =10,
     &          nbli =11,nmwt =12,nmwp =13,nmwu =14,nmwv =15,
     &          nzs  =16,nslp =17,nrh1 =18,nrh2 =19,nrh3 =20,
     &          ns1th=21,ns1t =22,ns1o =23,ns1r =24,ns1u =25,
     &          ns1v =26,ns1q =27,ntcl =28,ntpr=29)
c os,o,r,qs
      dimension kslp(2)
      logical lpcl,lppr,lpo3
      dimension si(im,km+1),sl(im,km)
      dimension zs(im),ps(im),os(im),clat(im)
      dimension u(ix,km),v(ix,km),o(im,km)
      dimension t(ix,km),r(im,km),q(ix,km),qs(im,km)
      dimension cl(ix,km),pr(ix,km)
      dimension tpt(im,kt),qpt(im,kt)
      dimension zm(im,*),tm(im,*)
      dimension sun(im,nsun)
      dimension wrk(im)
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 ,cp= 1.0046e+3 )
      parameter(rocp=rd/cp)
      parameter(pm1=1.e5,tm1=287.45,zm1=113.,zm2=5572.)
      parameter(fslp=g*(zm2-zm1)/(rd*tm1))
      parameter(strh1=0.44,strh2=0.72,strh3=0.44,
     &          sbrh1=1.00,sbrh2=0.94,sbrh3=0.72)
      parameter(sl1=0.9950)
c
c  surface orography, surface pressure and surface pressure tendency
c
      do i=1,im
        sun(i,nzs)=zs(i)
        sun(i,nps)=ps(i)*1.e3
        sun(i,npst)=os(i)
      enddo
c 
c  column precipitable water and relative humidity
c
      do i=1,im
        sun(i,ntpw)=0.
        wrk(i)=0.
      enddo
      do k=1,km
        do i=1,im
          ds=si(i,k)-si(i,k+1)
          sun(i,ntpw)=sun(i,ntpw)+q(i,k)*ds
          wrk(i)=wrk(i)+qs(i,k)*ds
        enddo
      enddo
      do i=1,im
        sun(i,ntrh)=min(max(sun(i,ntpw)/wrk(i),0.),1.)*100.
        sun(i,ntpw)=sun(i,ntpw)*ps(i)*1.e3/g
      enddo
c
c  sundry tropopause fields
c
      call sig2tp(im,ix,km,clat,sl,ps,u,v,t,
     &            sun(1,ntpp),sun(1,ntpu),sun(1,ntpv),
     &            sun(1,ntpt),sun(1,ntpsh))
      do i=1,im
        sun(i,ntpp)=sun(i,ntpp)*1.e3
      enddo
c 
c  lifted index
c
      if(kli.gt.0) then
        call liftix(im,ix,kt,pt,ps,tpt,qpt,tm(1,kli),
     &              sun(1,nsli),sun(1,nbli))
      else
        do i=1,im
          sun(i,nsli)=0.
          sun(i,nbli)=0.
        enddo
      endif
c
c  sundry maxwind fields
c
      call sig2mw(im,ix,km,sl,ps,u,v,t,
     &            sun(1,nmwp),sun(1,nmwu),sun(1,nmwv),sun(1,nmwt))
      do i=1,im
        sun(i,nmwp)=sun(i,nmwp)*1.e3
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  sea level pressure
      if(kslp(1).gt.0.and.kslp(2).gt.0) then
        k1=kslp(1)
        k2=kslp(2)
        do i=1,im
          sun(i,nslp)=pm1*exp(fslp*zm(i,k1)/(zm(i,k2)-zm(i,k1)))
        enddo
      else
        do i=1,im
          sun(i,nslp)=0.
        enddo
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  average relative humidity 1
      do i=1,im
        sun(i,nrh1)=0.
        wrk(i)=0.
      enddo
      do k=1,km
        do i=1,im
          ds=min(si(i,k),sbrh1)-max(si(i,k+1),strh1)
          if(ds.gt.0.) then
            sun(i,nrh1)=sun(i,nrh1)+q(i,k)*ds
            wrk(i)=wrk(i)+qs(i,k)*ds
          endif
        enddo
      enddo
      do i=1,im
        sun(i,nrh1)=min(max(sun(i,nrh1)/wrk(i),0.),1.)*100.
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  average relative humidity 2
      do i=1,im
        sun(i,nrh2)=0.
        wrk(i)=0.
      enddo
      do k=1,km
        do i=1,im
          ds=min(si(i,k),sbrh2)-max(si(i,k+1),strh2)
          if(ds.gt.0.) then
            sun(i,nrh2)=sun(i,nrh2)+q(i,k)*ds
            wrk(i)=wrk(i)+qs(i,k)*ds
          endif
        enddo
      enddo
      do i=1,im
        sun(i,nrh2)=min(max(sun(i,nrh2)/wrk(i),0.),1.)*100.
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  average relative humidity 3
      do i=1,im
        sun(i,nrh3)=0.
        wrk(i)=0.
      enddo
      do k=1,km
        do i=1,im
          ds=min(si(i,k),sbrh3)-max(si(i,k+1),strh3)
          if(ds.gt.0.) then
            sun(i,nrh3)=sun(i,nrh3)+q(i,k)*ds
            wrk(i)=wrk(i)+qs(i,k)*ds
          endif
        enddo
      enddo
      do i=1,im
        sun(i,nrh3)=min(max(sun(i,nrh3)/wrk(i),0.),1.)*100.
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  bottom sigma fields interpolated from first two model layers
      do i=1,im
        f2=log(sl(i,1)/sl1)/log(sl(i,1)/sl(i,2))
        sl1k=(sl1*1.e-2)**(-rocp)
        sun(i,ns1t)=t(i,1)+f2*(t(i,2)-t(i,1))
        sun(i,ns1th)=sun(i,ns1t)*sl1k*ps(i)**(-rocp)
        sun(i,ns1o)=o(i,1)+f2*(o(i,2)-o(i,1))
        sun(i,ns1r)=r(i,1)+f2*(r(i,2)-r(i,1))
        sun(i,ns1u)=u(i,1)+f2*(u(i,2)-u(i,1))
        sun(i,ns1v)=v(i,1)+f2*(v(i,2)-v(i,1))
        sun(i,ns1q)=q(i,1)+f2*(q(i,2)-q(i,1))
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  total cloud water
      if(lpcl) then
        do i=1,im
          sun(i,ntcl)=0.
        enddo
        do k=1,km
          do i=1,im
            ds=si(i,k)-si(i,k+1)
            sun(i,ntcl)=sun(i,ntcl)+cl(i,k)*ds
          enddo
        enddo
        do i=1,im
          sun(i,ntcl)=max(sun(i,ntcl),0.)*ps(i)*1.e3/g
        enddo
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  total precipatation water
      if(lppr) then
        do i=1,im
          sun(i,ntpr)=0.
        enddo
        do k=1,km
          do i=1,im
            ds=si(i,k)-si(i,k+1)
            sun(i,ntpr)=sun(i,ntpr)+pr(i,k)*ds
          enddo
        enddo
        do i=1,im
          sun(i,ntpr)=max(sun(i,ntpr),0.)*ps(i)*1.e3/g
        enddo
      endif
      return
      end
