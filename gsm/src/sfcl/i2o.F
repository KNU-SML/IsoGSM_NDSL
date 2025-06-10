      subroutine i2o(io,j,proj ,orient ,truth ,cotru ,
     1               delx ,dely ,rlat1 ,rlon1 ,ii, jj,
     2               projo,oriento,trutho,cotruo,
     3               delxo,delyo,rlat1o,rlon1o,glao,
     4               x00,y00,x00o,y00o,
     5               n00,n10,n11,n01,d00,d10,d11,d01,
     6               glon,glat,gaul)
#include <paramodel.h>
c
      real    d00(io),d10(io),d11(io),d01(io)
      integer n00(io),n10(io),n11(io),n01(io)
c
      real glat(jj),glon(ii),gaul(jj)
c
      logical loutside
c
      real, allocatable :: rlato(:),rlono(:)
      real, allocatable :: xlat (:), xlon(:)
      real, allocatable :: xo   (:),   yo(:)
      real, allocatable :: x    (:),    y(:)
c
      allocate (rlato(io),rlono(io))
      allocate (xlat (io),xlon (io))
      allocate (x    (io),y    (io))
c
c  get output latitude
c
      if(abs(projo).eq.3..or.abs(projo).eq.4.) then
        if(projo.eq.3.) then
          do i=1,io
            rlato(i)=rlat1o+delyo*float(j-1)/1000.
c           rlato(i)=rlat1o+delyo*float(j-1)
          enddo
        elseif(abs(projo).eq.4.) then
          do i=1,io
            rlato(i)=glao
          enddo
        elseif(projo.eq.-3.) then
          do i=1,io
            rlato(i)=rlat1o-delyo*float(j-1)/1000.
c           rlato(i)=rlat1o-delyo*float(j-1)
          enddo
        endif
c
c  get output longitude (and latitude)
c
        do i=1,io
          rlono(i)=mod(rlon1o+delxo*float(i-1)/1000.+360.,360.)
c         rlono(i)=mod(rlon1o+delxo*float(i-1)+360.,360.)
        enddo
        if(abs(proj).eq.3..or.abs(proj).eq.4.) then
          call ll2gg(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1,
     1               rlato,rlono,x,y)
        else
          call ll2xyr(io,proj,orient,truth,cotru,
     1                rlato,rlono,x,y)
        endif
c
      else
        allocate (xo(io),yo(io))
        do i=1,io
          yo(i)=y00o+(j-1)*delyo
        enddo
        do i=1,io
          xo(i)=x00o+(i-1)*delxo
        enddo
        call xy2llr(io,projo,oriento,trutho,cotruo,xo,yo,
     1             rlato,rlono)
        if(abs(proj).eq.3..or.abs(proj).eq.4.) then
          call ll2gg(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1,
     1               rlato,rlono,x,y)
        else
          call ll2xyr(io,proj,orient,truth,cotru,rlato,rlono,
     1                x,y)
        endif
        deallocate (xo,yo)
      endif
c
c  input latitude and longitude
c
      if(abs(proj).eq.3..or.abs(proj).eq.4.) then
        do i=1,io
          xlon(i)=x(i)
          xlat(i)=y(i)
        enddo
      else
        do i=1,io
          xlon(i)=(x(i)-x00)/delx+1
          xlat(i)=(y(i)-y00)/dely+1
        enddo
      endif
c
c  computation of coefficients
c
      do i=1,io
        loutside=.false.
c       if(xlon(i).lt.1) loutside=.true.
c       if(xlat(i).lt.1) loutside=.true.
c       if(xlon(i).gt.ii) loutside=.true.
c       if(xlat(i).gt.jj) loutside=.true.
        if(xlon(i).lt.0) loutside=.true.
        if(xlat(i).lt.0) loutside=.true.
        if(xlon(i).gt.ii+1) loutside=.true.
        if(xlat(i).gt.jj+1) loutside=.true.
        lon=max(xlon(i),1.)
        lat=max(xlat(i),1.)
        lon=min(lon,ii)
        lat=min(lat,jj)
        d00(i)=max(min((1.-(xlon(i)-lon)) * (1.-(xlat(i)-lat)),1.),0.)
        d10(i)=max(min((xlon(i)-lon) * (1.-(xlat(i)-lat)),1.),0.)
        d11(i)=max(min((xlon(i)-lon) * (xlat(i)-lat),1.),0.)
        d01(i)=max(min((1.-(xlon(i)-lon)) * (xlat(i)-lat),1.),0.)
        n00(i)=lon   + (lat   -1)*ii
        n10(i)=lon+1 + (lat   -1)*ii
        n11(i)=lon+1 + (lat+1 -1)*ii
        n01(i)=lon   + (lat+1 -1)*ii
        if(loutside) then
          d00(i)=-999.
          d10(i)=-999.
          d11(i)=-999.
          d01(i)=-999.
          n00(i)=-999
          n10(i)=-999
          n11(i)=-999
          n01(i)=-999
        endif
      enddo
c
      deallocate (rlato,rlono,xlat,xlon,x,y)
      return
      end
