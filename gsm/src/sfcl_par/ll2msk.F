      subroutine ll2msk(lugb,fnmskg,
     &                  inttyp,rout,io,jo,
     &                  projo,oriento,trutho,cotruo,
     &                  delxo,delyo,rlat1o,rlat2o,rlon1o,rlon2o)
#include <paramodel.h>
#undef VECTOR
#ifdef NEC
#define VECTOR
#endif
#ifdef ES
#define VECTOR
#endif
cfpp$ noconcur r
c
c  interpolation from lat/lon or gaussian grid to other lat/lon grid
c
c  inttyp ... type of interpolation.  0 = bilinear
c                                     1 = take the closest point value
c
      dimension rout(io*jo)
c
      parameter(imxgrmsk_max=2280,jmxgrmsk_max=1141)
c
      character*128 fnmskg
      integer kpds_grmsk(25)
      data kpds_grmsk/4*-1, 81, -1,  -1,18*-1/
c
      real, allocatable :: grmsk(:)
      real, allocatable :: glon(:),glat(:),gaul(:)
      real, allocatable :: glao(:)
c
      real,    allocatable ::  d00(:),d10(:),d11(:),d01(:)
      integer, allocatable ::  n00(:),n10(:),n11(:),n01(:)
      integer, allocatable ::  ij1(:)
c
      allocate (d00(io),d10(io),d11(io),d01(io))
      allocate (n00(io),n10(io),n11(io),n01(io))
      allocate (ij1(io))
c
c  read in high resolution global lat/lon land sea mask
c
      allocate (grmsk(imxgrmsk_max*jmxgrmsk_max))
      print *,'----------- start reading hires lat/lon mask -----'
      call fixrd(lugb,fnmskg,kpds_grmsk,
     &           0,0,0,0,0.,
     &           imxgrmsk_max,jmxgrmsk_max,
     &           imxgrmsk,jmxgrmsk,
     &           kkpds5,kkgds1,
     &           proji,orienti,truthi,cotrui,
     &           delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i,
     &           dlon,dlat,wlon,rnlat,
     &           grmsk,iret)
      if(iret.gt.0) then
        print *,'high res global lat/lon land sea mask read failed.'
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
      if(imxgrmsk_max*jmxgrmsk_max.lt.imxgrmsk*jmxgrmsk) then
        print *,'dimension size overflow in ll2msk'
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
      print *,'----------- finish reading hires lat/lon mask -----'
      close(lugb)
c
#ifdef DBG
      print *,' --- into ll2msk --- '
      print *,'ll2msk:grmsk'
      call nntprt(grmsk,imxgrmsk,jmxgrmsk)
#endif
c
      ijo = io*jo
      if(imxgrmsk.eq.1.or.jmxgrmsk.eq.1) then
        print*,' wrong imxgrmsk or jmxgrmsk in ll2msk',imxgrmsk,jmxgrmsk
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
      if(delxi.eq.0..or.delyi.eq.0.) then
        print *,'delxi or delyi is zero in ll2msk.'
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
#ifdef DBG
      print *,'delxi=',delxi
      print *,'delyi=',delyi
      print *,'rlon1i=',rlon1i
      print *,'rlat1i=',rlat1i
      print *,' prepare grid to grid interpolation.'
      print *,' input grid : '
      print *,' proj=',proji,' orient=',orienti,
     1        ' truth=',truthi,' cotru=',cotrui,
     1        ' delx=',delxi,' dely=',delyi,
     1        ' rlat1=',rlat1i,' rlon1=',rlon1i
      print *,' output grid : '
      print *,' projo=',projo,' oriento=',oriento,
     1        ' trutho=',trutho,' cotruo=',cotruo,
     1        ' delxo=',delxo,' delyo=',delyo,
     1        ' rlat1o=',rlat1o,' rlon1o=',rlon1o
#endif
c
      allocate (glon(imxgrmsk),glat(jmxgrmsk),gaul(jmxgrmsk))
      allocate (glao(jo))
c
      call i2oini(proji,orienti,truthi,cotrui,
     1            delxi,delyi,rlat1i,rlon1i,imxgrmsk,jmxgrmsk,
     2            projo,oriento,trutho,cotruo,
     3            x00,y00,x00o,y00o,
     4            rlat1o,rlon1o,jo,glao,glon,glat,gaul)
c
c  inttyp=1  take the closest point value
c
#ifdef DBG
      print *,'in ll2msk.  inttyp=',inttyp
#endif
c
      do j=1,jo
        call i2o(io,j,proji,orienti,truthi,cotrui,
     1           delxi,delyi,rlat1i,rlon1i,imxgrmsk,jmxgrmsk,
     2           projo,oriento,trutho,cotruo,
     3           delxo,delyo,rlat1o,rlon1o,glao(j),
     4           x00,y00,x00o,y00o,
     5           n00,n10,n11,n01,d00,d10,d11,d01,
     6           glon,glat,gaul)
c
        if(inttyp.eq.1) then
          do i=1,io
            weimax = 0.
            if(weimax.le.d00(i)) then
              weimax = d00(i)
              ij1(i) = n00(i)
            elseif(weimax.le.d01(i)) then
              weimax = d01(i)
              ij1(i) = n01(i)
            elseif(weimax.le.d10(i)) then
              weimax = d10(i)
              ij1(i) = n10(i)
            elseif(weimax.le.d11(i)) then
              weimax = d11(i)
              ij1(i) = n11(i)
            endif
          enddo
          do i=1,io
            ij=(j-1)*io+i
            rout(ij)=grmsk(ij1(i))
          enddo
        else
          do i=1,io
            ij=(j-1)*io+i
            a00 = grmsk(n00(i))
            a10 = grmsk(n10(i))
            a11 = grmsk(n11(i))
            a01 = grmsk(n01(i))
            wi1j1 = d00(i)
            wi2j1 = d10(i)
            wi1j2 = d01(i)
            wi2j2 = d11(i)
            wsum  = wi1j1+wi2j1+wi1j2+wi2j2
            wsumiv = 1./wsum
            rout(ij)=(wi1j1*a00+wi2j1*a10+
     1                wi1j2*a01+wi2j2*a11)*wsumiv
          enddo
        endif
      enddo
c
#ifdef DBG
      call rmaxmin(rout,io*jo,1,1,1,'mask for input grib from grmsk')
#endif
c
      deallocate (grmsk,glon,glat,gaul,glao)
      deallocate (d00,d10,d11,d01,n00,n10,n11,n01,ij1)
      return
      end
