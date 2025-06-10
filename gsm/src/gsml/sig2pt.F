      subroutine sig2pt(im,ix,km,si,sl,
     &                  ps,us,vs,ts,qs,qss,
     &                  kt,pt,upt,vpt,tpt,qpt,rpt)
c$$$  subprogram documentation block
c
c  subprogram:    sig2pt      sigma to pressure thickness
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: interpolates winds, temperature and humidity
c   from the sigma coordinate system to constant pressure thicknesses
c   above the ground.
c
c program history log:
c   94-07-08  iredell
c
c usage:    call sig2pt(im,ix,km,si,sl,
c    &                  ps,us,vs,ts,qs,qss,
c    &                  kt,pt,upt,vpt,tpt,qpt,rpt)
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
c     ts       - real (ix,km) temperature in k
c     qs       - real (ix,km) specific humidity in kg/kg
c     qss      - real (im,km) saturated specific humidity in kg/kg
c     kt       - integer number of pressure thickness layers
c     pt       - real pressure thickness in kpa
c
c   output argument list:
c     upt      - real (im,kt) zonal wind in m/s
c     vpt      - real (im,kt) merid wind in m/s
c     tpt      - real (im,kt) temperature in k
c     qpt      - real (im,kt) specific humidity in kg/kg
c     rpt      - real (im,kt) relative humidity in percent
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension si(im,km+1),sl(im,km),ps(im)
      dimension us(ix,km),vs(ix,km),ts(ix,km),qs(ix,km),qss(im,km)
      dimension upt(im,kt),vpt(im,kt),tpt(im,kt),qpt(im,kt),rpt(im,kt)
c
      rmin=1.e+33
      do i=1,im
        rmin=min(rmin,ps(i))
      enddo
      rmin=rmin/pt
      rmax=-1.e+33
      do i=1,im
        rmax=max(rmax,ps(i))
      enddo
      rmax=rmax/pt
c
      do k=1,kt
        do i=1,im
          upt(i,k)=0.
          vpt(i,k)=0.
          tpt(i,k)=0.
          qpt(i,k)=0.
          rpt(i,k)=0.
        enddo
        do ks=1,km
          do i=1,im
            if(k-(1-si(i,ks))*rmin.gt.0.
     $            .and.k-(1-si(i,ks+1))*rmax.lt.1.) then
              r=ps(i)/pt
              fks=min(k-(1-si(i,ks))*r,1.)-max(k-(1-si(i,ks+1))*r,0.)
              if(fks.gt.0.) then
                upt(i,k)=upt(i,k)+fks*us(i,ks)
                vpt(i,k)=vpt(i,k)+fks*vs(i,ks)
                tpt(i,k)=tpt(i,k)+fks*ts(i,ks)
                qpt(i,k)=qpt(i,k)+fks*qs(i,ks)
                rpt(i,k)=rpt(i,k)+fks*qss(i,ks)
              endif
            endif
          enddo
        enddo
        do i=1,im
          rpt(i,k)=min(max(qpt(i,k)/rpt(i,k),0.),1.)*100.
        enddo
      enddo
c
      return
      end
