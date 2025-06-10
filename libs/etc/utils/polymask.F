      program polymask
c
      parameter(nptsmax=100)
      real x(nptsmax),y(nptsmax)
c
      real, allocatable :: mask(:,:),masko(:,:)
      real*4, allocatable :: mask4(:,:)
c
      data x/nptsmax*-999./
      data y/nptsmax*-999./
c
      parameter(idim=500,jdim=250)
c    
      namelist/nampoly/ x,y,
     1                  io,jo,proj,truth,orient,delx,dely,
     2                  cenlat,cenlon,rftgrd,btmgrd 
c
      read(5,nampoly)
c
      do n=1,nptsmax 
        if(x(n).eq.-999.) then
         npts=n-1
         go to 1
        endif
      enddo
    1 continue
      print *,'npts=',npts
c
      do n=1,npts
        if(x(n).lt.0.) x(n)=x(n)+360.
      enddo
c
      xmax=-180.
      xmin=360.
      do n=1,npts
        xmax=max(xmax,x(n))
        xmin=min(xmin,x(n))
      enddo
      do n=1,npts
        if(xmax.eq.x(n)) imax=n
        if(xmin.eq.x(n)) imin=n
      enddo
c
      ymax=-90.
      ymin=90.
      do n=1,npts
        ymax=max(ymax,y(n))
        ymin=min(ymin,y(n))
      enddo
      do n=1,npts
        if(ymax.eq.y(n)) jmax=n
        if(ymin.eq.y(n)) jmin=n
      enddo
c
      safefac=5.
      xsafe=delx/1000./111.*safefac
      ysafe=dely/1000./111.*safefac
      glonw=max(xmin-xsafe,0.)
      glats=max(ymin-ysafe,0.)
      glone=min(xmax+xsafe,360.)
      glatn=min(ymax+ysafe,90.)
      dlon=(glone-glonw)/idim
      dlat=(glatn-glats)/(jdim-1)
c
      if(npts.gt.nptsmax) then
        print *,'npts .gt.',nptsmax
        call abort
      endif
c
      allocate (mask(idim,jdim))
      allocate (masko(io,jo))
      allocate (mask4(io,jo))
c
      do j=1,jdim
        do i=1,idim
          mask(i,j)=0
        enddo
      enddo
c
c  compute the longitude where the polygon crosses the i-grid line
c
      do j=1,jdim
        flat=glats+dlat*float(j-1)
        if(flat.lt.ymin.or.flat.gt.ymax) go to 10
        i1=0
        i2=0
        do n=1,npts
          nn=jmin+n-1
          nnp=nn+1
          if(nn.gt.npts) nn=nn-npts
          if(nnp.gt.npts) nnp=nnp-npts
          if(i1.eq.0.and.flat.ge.y(nn).and.flat.lt.y(nnp)) then
             xcr=x(nn)+(x(nnp)-x(nn))/(y(nnp)-y(nn))*(flat-y(nn))
             i1=(xcr-glonw)/dlon+1
          endif
          nn=jmin-n+1
          nnm=nn-1
          if(nn .lt.1) nn =npts+nn
          if(nnm.lt.1) nnm=npts+nnm
          if(i2.eq.0.and.flat.ge.y(nn).and.flat.lt.y(nnm)) then
             xcr=x(nn)+(x(nnm)-x(nn))/(y(nnm)-y(nn))*(flat-y(nn))
             i2=(xcr-glonw)/dlon+1
          endif
        enddo
        if(i1.ne.0.and.i2.ne.0) then
          do i=min(i1,i2),max(i1,i2)
            mask(i,j)=1
          enddo
        endif
   10 continue
      enddo
c
c  compute the longitude where the polygon crosses the i-grid line
c
      do i=1,idim
        flon=glonw+dlon*float(i-1)
        if(flon.lt.xmin.or.xlon.gt.xmax) go to 20
        j1=0
        j2=0
        do n=1,npts
          nn=imin+n-1
          nnp=nn+1
          if(nn.gt.npts) nn=nn-npts
          if(nnp.gt.npts) nnp=nnp-npts
          if(j1.eq.0.and.flon.ge.x(nn).and.flon.lt.x(nnp)) then
             ycr=y(nn)+(y(nnp)-y(nn))/(x(nnp)-x(nn))*(flon-x(nn))
             j1=(ycr-glats)/dlat+1
          endif
          nn=imin-n+1
          nnm=nn-1
          if(nn .lt.1) nn =npts+nn
          if(nnm.lt.1) nnm=npts+nnm
          if(j2.eq.0.and.flon.ge.x(nn).and.flon.lt.x(nnm)) then
             ycr=y(nn)+(y(nnm)-y(nn))/(x(nnm)-x(nn))*(flon-x(nn))
             j2=(ycr-glats)/dlat+1
          endif
        enddo
        if(j1.ne.0.and.j2.ne.0) then
          do j=min(j1,j2),max(j1,j2)
            mask(i,j)=mask(i,j)+1
          enddo
        endif
   20   continue
      enddo
c
      do j=1,jdim
        do i=1,idim
          mask(i,j)=max(mask(i,j)-1.,0.)
        enddo
      enddo
c
c     call nntprt(mask,idim,jdim)
c
      call ll2rsm(mask,idim,jdim,glonw,glats,dlon,dlat,
     1            proj,truth,orient,delx,dely,
     2            cenlat,cenlon,rftgrd,btmgrd,
     3            masko,io,jo)
c
c     call nntprt(masko,io,jo)
c
      do j=1,jo
        do i=1,io
          mask4(i,j)=masko(i,j)
        enddo
      enddo
      open(1,file='polymask.data',form='unformatted')
      write(1) mask4
      call makectls(io,jo,
     1              proj,truth,orient,delx,dely,
     2              cenlat,cenlon,rftgrd,btmgrd,2)
c
      stop
      end
      subroutine nntprt(data,imax,jmax)
      dimension data(imax*jmax)
c
      rmax=data(1)
      rmin=data(1)
      do ij=1,imax*jmax
        rmax=max(rmax,data(ij))
        rmin=min(rmin,data(ij))
      enddo
      if(rmax.eq.rmin) then
        print *,'constant field. value=',rmax
        return
      endif
      rmxmn=rmax-rmin
      if(rmax.ne.0.) then
        rn=log10(abs(rmax))
        if(rn.lt.0.) then
          n=rn-1.0
        else
          n=rn
        endif
      else
        n=-9999
      endif
      if(rmin.ne.0.) then
        rm=log10(abs(rmin))
        if(rm.lt.0.) then
          m=rm-1.0
        else
          m=rm
        endif
      else
        m=-9999
      endif
      n=max(n,m)
      fact=10.**(-n)
      print *,'nntprt:rmax=',rmax,' rmin=',rmin,' fact=',fact
c
      ilast=0
      i1=1
      i2=80
      if(i2.ge.imax) then
        ilast=1
        i2=imax
      endif
 1112 continue
      write(6,*) ' '
c     do j=1,jmax
      do j=jmax,1,-1
        write(6,1111) (nint(data(imax*(j-1)+i)*fact),i=i1,i2)
      enddo
      if(ilast.eq.1) return
      i1=i1+80
      i2=i1+79
      if(i2.ge.imax) then
        ilast=1
        i2=imax
      endif
      go to 1112
 1111 format(80i1)
      return
      end
      subroutine ll2rsm(rdat,imin,jmin,glonw,glats,dlon,dlat,
     1            proj,truth,orient,delx,dely,
     2            cenlat,cenlon,rftgrd,btmgrd,
     2            rout,io,jo)
cfpp$ noconcur r
c
c  interpolation from global lat/lon grid to other grid
c
      parameter(pi=3.1415926536)
      dimension rdat(imin*jmin)
      dimension rout(io*jo)
c
      real,    allocatable ::  d00(:),d10(:),d11(:),d01(:)
      integer, allocatable ::  n00(:),n10(:),n11(:),n01(:)
c
      allocate (d00(io),d10(io),d11(io),d01(io))
      allocate (n00(io),n10(io),n11(io),n01(io))
c
c  find lat/lon of the bottom left corner of the domain
c
      call setgrd(rlat1,rlon1,
     1            proj,truth,orient,delx,dely,
     2            cenlat,cenlon,rftgrd,btmgrd,dlamda0)
c
      print *, ' prepare global lat/lon to output grid interpolation.'
      print *, ' output grid : '
      print *, ' proj=',proj,' orient=',orient,
     1         ' truth=',truth,' cenlat=',cenlat,' cenlon=',cenlon,
     2         ' delx=',delx,' dely=',dely,
     3         ' rftgrd=',rftgrd,' btmgrd=',btmgrd,
     5         ' rlat1=',rlat1,' rlon1=',rlon1
c
      do j=1,jo
        call i2o(io,j,
     1           glonw,glats,dlon,dlat,imin,jmin,
     2           proj,orient,truth,cotru,
     3           delx,dely,rlat1,rlon1,
     5           n00,n10,n11,n01,d00,d10,d11,d01)
c
c  bilinear interpolation without mask
c
        do i=1,io
          ij=(j-1)*io+i
          a00 = rdat(n00(i))
          a10 = rdat(n10(i))
          a11 = rdat(n11(i))
          a01 = rdat(n01(i))
          wi1j1 = d00(i)
          wi2j1 = d10(i)
          wi1j2 = d01(i)
          wi2j2 = d11(i)
          wsum  = wi1j1+wi2j1+wi1j2+wi2j2
          wsumiv = 1./wsum
          rout(ij)=(wi1j1*a00+wi2j1*a10+
     1              wi1j2*a01+wi2j2*a11)*wsumiv
        enddo
      enddo
c
      deallocate (d00,d10,d11,d01,n00,n10,n11,n01)
      return
      end
      subroutine setgrd(flat1,flon1,
     1                  rproj,rtruth,rorient,rdelx,rdely,
     2                  rcenlat,rcenlon,rlftgrd,rbtmgrd,dlamda0)
c
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  setgrd
c   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
c   prgmmr:  yifeng cui, mp print       org: sdsc       date: 02-09-11
c
c abstract:  setup grid latitude, longitude, and map factor etc for
c            regional grid.
c
c program history log:
c
c usage:    call setgrd(flat1,flon1,delx,dely,dlamda0)
c   output argument list:
c     flat	- latitude of left bottom corner
c     flon	- logitude of left bottom corner
c     delx	- grid spacing in x direction
c     dely	- grid spacing in y direction
c     dlamda0	- angle between 0 longitude and x positive axies
c
c   input argument list:
c
c   input files: none
c
c   output files: none
c
c   subprograms called: none
c
c   remark: none
c
c attributes:
c   language: fortran 77.
c   machine:  cray ymp.
c
c$$$
      parameter(pi=3.141593,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.)
      parameter(rerth=6.3712e+6)
c
      delx = rdelx
      dely = rdely
c
c --------- setup regional lat/lon and map factor -----
c
c if proj=0  do mercater projection
c if proj=1  do north polar projection
c if proj=-1 do south polar projection
c
      nproj = rproj
c
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
c ++++++++++++++++++++++++++++++++++++++
c polar projection
c ++++++++++++++++++++++++++++++++++++++
      truth  = rtruth * rad
      truth  = nproj * truth
      orient  = rorient * rad
      dlamda0 = orient + hfpi
      a2 =  rerth * ( 1.0 + sin(truth) )
      radlat = rcenlat * rad
      radlon = rcenlon * rad - dlamda0
      radlat = nproj * radlat
      radlon = nproj * radlon
      yyy = a2 * cos( radlat )/(1. + sin( radlat ) )
      cenlat = rcenlat
      if( abs(cenlat) .eq. 90. ) then yyy = 0.0
      y00 = yyy * sin( radlon ) - ( rbtmgrd -1.) * dely
      x00 = yyy * cos( radlon ) - ( rlftgrd -1.) * delx
c
      j=1
      ys = y00 + (j-1)*dely
      i=1
      x = x00 + (i-1)*delx
      if( x .gt. 0. e 0 ) then
         flons = atan(ys/x)
      else if ( x .lt. 0. e 0 ) then
         flons = pi + atan(ys/x)
      else
         flons = hfpi
         if( ys .lt. 0. e 0 ) flons = flons * 3.0
      endif
      flons = nproj * flons + dlamda0
      flons = mod(flons,twopi)
      if(flons.lt.0. e 0) flons = twopi + flons
c
      rsoa2 = sqrt( x*x + ys*ys )/a2
      flats = hfpi - 2.0 e 0 * atan(rsoa2)
      flat1 = nproj * flats
      flon1 = flons
c
      else if ( nproj.eq.0 ) then
c
c ++++++++++++++++++++++++++++
c mercator
c ++++++++++++++++++++++++++++
      truth  = rtruth * rad
      cenlat = rcenlat * rad
      cenlon = rcenlon * rad
      cenlon = mod(cenlon,twopi)
      if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
      a2 =  rerth * cos( truth )
      x0 = 0.0
      y0 = a2 * log( abs( tan( qtpi + 0.5 * cenlat ) ) )
      x00 = -( rlftgrd - 1. ) * delx - 0.5 * delx
      y00 = -( rbtmgrd - 1. ) * dely - 0.5 * dely
      dlamda0 = 0.0
c
      j=1
      lats = j
      ys = y00 + (j-1)*dely + y0
c
      i=1
      x = x00 + (i-1)*delx + x0
      flons = x / a2 + cenlon
      flons = mod(flons,twopi)
      if(flons.lt.0. e 0) flons = twopi + flons
c
      flats = 2.0 *( atan( exp( ys/a2 ) ) - qtpi )
      flat1 = flats*180./pi
      flon1 = flons*180./pi
c
      endif
c
      return
      end
      subroutine i2o(io,j,
     1               glonw,glats,dlon ,dlat ,ii, jj,
     2               projo,oriento,trutho,cotruo,
     3               delxo,delyo,rlat1o,rlon1o,
     5               n00,n10,n11,n01,d00,d10,d11,d01)
c
      real    d00(io),d10(io),d11(io),d01(io)
      integer n00(io),n10(io),n11(io),n01(io)
c
      real, allocatable :: glao(:)
      real, allocatable :: glat(:),glon(:)
      real, allocatable :: rlato(:),rlono(:)
      real, allocatable :: xlat (:), xlon(:)
      real, allocatable :: xo   (:),   yo(:)
      real, allocatable :: x    (:),    y(:)
c
      data ifp/0/
      save ifp
c
      allocate (glat (jj),glon (ii))
      allocate (rlato(io),rlono(io))
      allocate (xlat (io),xlon (io))
      allocate (x    (io),y    (io))
c
c  get input lat/lon
c
      if(ifp.eq.0) then 
        do jy = 1,jj
          glat(jy) = glats + dlat * float(jy-1)
        enddo
        do i = 1,ii
          glon(i) = glonw+dlon * float(i-1)
          if(glon(i).lt.0.) glon(i) = glon(i) + 360.
        enddo
        ifp=1
      endif
c
c  get output latitude
c
      if(abs(projo).eq.3..or.abs(projo).eq.4.) then
        if(projo.eq.3.) then
          do i=1,io
            rlato(i)=rlat1o+delyo*float(j-1)/1000.
          enddo
        elseif(abs(projo).eq.4.) then
          allocate (glao(jj))
          call gaulat(glao,jj)
          do jy = 1,jj
            if(projo.eq.-4.) then  
               glao(jy) = 90.-glao(jy)
            elseif(projo.eq.4) then
               glao(jj-jy+1) = 90.-glao(jy)
            endif
          enddo
          do i=1,io
            rlato(i)=glao(j)
          enddo
          deallocate (glao)
        elseif(projo.eq.-3.) then
          do i=1,io
            rlato(i)=rlat1o-delyo*float(j-1)/1000.
          enddo
        endif
c
c  get output longitude
c
        do i=1,io
          rlono(i)=mod(rlon1o+delxo*float(i-1)/1000.+360.,360.)
        enddo
c
      else
        call ll2xyr(projo,oriento,trutho,cotruo,rlat1o,rlon1o,
     1              x00o,y00o)
        allocate (xo(io),yo(io))
        do i=1,io
          yo(i)=y00o+(j-1)*delyo
        enddo
        do i=1,io
          xo(i)=x00o+(i-1)*delxo
        enddo
        call xy2llr(io,projo,oriento,trutho,cotruo,xo,yo,rlato,rlono)
        deallocate (xo,yo)
      endif
c
      call ll2gg(io,3.,glat,glon,ii,jj,dlon,dlat,rlato,rlono,xlon,xlat)
c
c  computation of coefficients
c
      do i=1,io
        lon=max(xlon(i),1.)
        lat=max(xlat(i),1.)
        lon=min(lon,ii)
        lat=min(lat,jj)
        d00(i)=(1.-(xlon(i)-lon)) * (1.-(xlat(i)-lat))
        d10(i)=(xlon(i)-lon) * (1.-(xlat(i)-lat))
        d11(i)=(xlon(i)-lon) * (xlat(i)-lat)
        d01(i)=(1.-(xlon(i)-lon)) * (xlat(i)-lat)
        n00(i)=lon   + (lat   -1)*ii
        n10(i)=lon+1 + (lat   -1)*ii
        n11(i)=lon+1 + (lat+1 -1)*ii
        n01(i)=lon   + (lat+1 -1)*ii
      enddo
c
      deallocate (glat,glon,rlato,rlono,xlat,xlon,x,y)
      return
      end
      subroutine xy2llr(io,cproj,corient,ctruth,ccotru,cx,cy,clat,clon)
c
      parameter(pi=3.1415926536,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.,rerth=6.3712e+6)
c
      real cx(io),cy(io)
      real clat(io),clon(io)
c
c input are all degree and output to global x y, not regional x y
c
c if proj=0   mercator projection
c if proj=1   north polar projection
c if proj=-1  south polar projection
c if proj=2   north lambert projection
c if proj=-2  south lambert projection
c
c polar projection
c
      nproj = cproj
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
        truth  = ctruth * rad
        truth  = nproj * truth
        orient  = corient * rad
        cenlon = mod(orient,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        dlamda0 = cenlon + hfpi
        a2 =  rerth * ( 1.0 + sin(truth) )
        radlat = 90. * rad
        radlon = 0.0 * rad - dlamda0
        radlat = nproj * radlat
        radlon = nproj * radlon
        do i=1,io
          x = cx(i)
          y = cy(i)
          if( x.gt.0.0 ) then
            blon = atan(y/x)
          else if( x.lt.0.0 ) then
            blon = pi + atan(y/x)
          else
            blon = hfpi
            if( y.lt.0.0 ) blon = blon * 3.0
          endif
          blon = blon + dlamda0
          blon = mod(blon,twopi)
          blon = nproj * blon
          rsoa2 = sqrt( x*x + y*y )/a2
          blat = hfpi - 2. * atan(rsoa2)
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
c
c  mercator
c
      else if ( nproj.eq.0 ) then
        truth  = ctruth * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        a2 =  rerth * cos( truth )
        dlamda0 = 0.0
        do i=1,io
          x = cx(i)
          y = cy(i)
          blon = x / a2 + cenlon
          blon = mod(blon,twopi)
          if(blon.lt.0. e 0) blon = twopi + blon
          blat = 2.*(atan(exp(y/a2))-qtpi)
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
c
c  lambert
c
      else if( nproj.eq.2 .or. nproj.eq.-2 ) then
        is=1
        if( nproj.lt.0 ) is=-1
        truth  = ctruth * rad
        cotru  = ccotru * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        if( ctruth.eq.ccotru ) then
          cone= cos (hfpi-is*truth)
        else
          cone=(log(cos(truth))-log(cos(cotru)))/
     1         (log(tan(qtpi-is*truth/2))-log(tan(qtpi-is*cotru/2)))
        endif
        dlamda0 = 0.0
        r00=rerth/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
        do i=1,io
          r = sqrt(cx(i)*cx(i) + cy(i)*cy(i) )
          blon = cenlon + asin(cx(i)/r) / cone
          blon = mod(blon+twopi,twopi)
          blat = hfpi - 2 * is * atan ( (r/r00)**(1./cone) )
          clat(i) = blat / rad
          clon(i) = blon / rad
        enddo
      endif
c
      return
      end
      subroutine gaulat(gaul,k)
c
      dimension a(2000)
      dimension gaul(1)
c
      esp=1.e-14
      c=(1.e0-(2.e0/3.14159265358979e0)**2)*0.25e0
      fk=k
      kk=k/2
      call bsslz1(a,kk)
      do 30 is=1,kk
      xz=cos(a(is)/sqrt((fk+0.5e0)**2+c))
      iter=0
   10 pkm2=1.e0
      pkm1=xz
      iter=iter+1
      if(iter.gt.10) go to 70
      do 20 n=2,k
      fn=n
      pk=((2.e0*fn-1.e0)*xz*pkm1-(fn-1.e0)*pkm2)/fn
      pkm2=pkm1
   20 pkm1=pk
      pkm1=pkm2
      pkmrk=(fk*(pkm1-xz*pk))/(1.e0-xz**2)
      sp=pk/pkmrk
      xz=xz-sp
      avsp=abs(sp)
      if(avsp.gt.esp) go to 10
      a(is)=xz
   30 continue
      if(k.eq.kk*2) go to 50
      a(kk+1)=0.e0
      pk=2.e0/fk**2
      do 40 n=2,k,2
      fn=n
   40 pk=pk*fn**2/(fn-1.e0)**2
   50 continue
      do 60 n=1,kk
      l=k+1-n
      a(l)=-a(n)
   60 continue
c
      radi=180.e0/(4.e0*atan(1.e0))
      do 211 n=1,k
      gaul(n)=acos(a(n))*radi
  211 continue
      return
   70 print *,'error in gauaw'
      stop
      end
      subroutine bsslz1(bes,n)
      save
c
      dimension bes(n)
      dimension bz(50)
c
      data pi/3.14159265358979e0/
      data bz         / 2.4048255577e0, 5.5200781103e0,
     $  8.6537279129e0,11.7915344391e0,14.9309177086e0,18.0710639679e0,
     $ 21.2116366299e0,24.3524715308e0,27.4934791320e0,30.6346064684e0,
     $ 33.7758202136e0,36.9170983537e0,40.0584257646e0,43.1997917132e0,
     $ 46.3411883717e0,49.4826098974e0,52.6240518411e0,55.7655107550e0,
     $ 58.9069839261e0,62.0484691902e0,65.1899648002e0,68.3314693299e0,
     $ 71.4729816036e0,74.6145006437e0,77.7560256304e0,80.8975558711e0,
     $ 84.0390907769e0,87.1806298436e0,90.3221726372e0,93.4637187819e0,
     $ 96.6052679510e0,99.7468198587e0,102.888374254e0,106.029930916e0,
     $ 109.171489649e0,112.313050280e0,115.454612653e0,118.596176630e0,
     $ 121.737742088e0,124.879308913e0,128.020877005e0,131.162446275e0,
     $ 134.304016638e0,137.445588020e0,140.587160352e0,143.728733573e0,
     $ 146.870307625e0,150.011882457e0,153.153458019e0,156.295034268e0/
c
      nn=n
      if(n.le.50) go to 12
      bes(50)=bz(50)
      do 5 j=51,n
    5 bes(j)=bes(j-1)+pi
      nn=49
   12 do 15 j=1,nn
   15 bes(j)=bz(j)
      return
      end
      subroutine ll2gg(io,proj,glat,glon,idim,jdim,delx,dely,
     1                 rlat,rlon,xlon,ylat)
c
      dimension glat(jdim),glon(idim)
      real rlat(io),rlon(io),xlon(io),ylat(io)
c
      do i=1,io
        do ix=2,idim
          if(glon(ix-1).le.rlon(i).and.glon(ix).gt.rlon(i)) then
             xlon(i) = ix-1 + (rlon(i)-glon(ix-1))/(glon(ix)-glon(ix-1))
          endif
        enddo
      enddo
c
      do i=1,io
        do j = 2,jdim
          if(glat(j-1).le.rlat(i).and.glat(j).gt.rlat(i)) then
            ylat(i) = j-1 + (rlat(i)-glat(j-1))/(glat(j)-glat(j-1))
          endif
        enddo
      enddo
c
      return
      end
      subroutine ll2xyr(cproj,corient,ctruth,ccotru,clat,clon,cx,cy)
c
      parameter(pi=3.1415926536,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.,rerth=6.3712e+6)
c
c input are all degree and output to global x y, not domain x y
c
c if proj=0  do mercater projection
c if proj=1  do north polar projection
c if proj=-1 do south polar projection
c if proj=2  do north lambert projection
c if proj=-2 do south lambert projection
c
      nproj = cproj
      blat = clat * rad
      blon = clon * rad
c
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
c ++++++++++++++++++++++++++++++++++++++
c polar projection
c ++++++++++++++++++++++++++++++++++++++
        truth  = ctruth * rad
        truth  = nproj * truth
        orient  = corient * rad
        cenlon = mod(orient,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        dlamda0 = cenlon + hfpi
        a2 =  rerth * ( 1.0 + sin(truth) )
        radlat = 90. * rad
        radlon = 0.0 * rad - dlamda0
        radlat = nproj * radlat
        radlon = nproj * radlon
c
c =========
        rsoa2 = tan( (hfpi-blat*nproj)*0.5 )
        x2py2 = ( rsoa2 * a2 ) ** 2.0
        blon = mod(blon,twopi)
        if(blon.lt.0. e 0) blon = twopi + blon
        rlon = nproj * (blon - dlamda0)
        rlon = amod(rlon,twopi)
        if( rlon.lt.0. ) rlon=twopi+rlon
        yox = tan(rlon)
        x = sqrt( x2py2/(1.+yox*yox) )
        y = sqrt( x2py2 - x*x )
        if( rlon.gt.hfpi .and. rlon.lt. pi+hfpi ) x = -x
        if( rlon.gt.pi .and. rlon.lt. twopi ) y = -y
c
      else if ( nproj.eq.0 ) then
c
c ++++++++++++++++++++++++++++
c do mercater
c ++++++++++++++++++++++++++++
        truth  = ctruth * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        a2 =  rerth * cos( truth )
        dlamda0 = 0.0
c
        blon = mod(blon,twopi)
        if(blon.lt.0. e 0) blon = twopi + blon
        x=a2*(blon-cenlon)
        y=a2*log(tan(blat/2.0+qtpi))
c
      else if( nproj.eq.2 .or. nproj.eq.-2 ) then
c
c ++++++++++++++++++++++++++++
c do lambert
c ++++++++++++++++++++++++++++
        is=1
        if( nproj.lt.0 ) is=-1
        truth  = ctruth * rad
        cotru  = ccotru * rad
        cenlon = corient * rad
        cenlon = mod(cenlon,twopi)
        if(cenlon.lt.0. e 0) cenlon = twopi + cenlon
        if( ctruth.eq.ccotru ) then
          cone= cos (hfpi-is*truth)
        else
          cone=(log(cos(truth))-log(cos(cotru)))/
     1         (log(tan(qtpi-is*truth/2))-log(tan(qtpi-is*cotru/2)))
        endif
        dlamda0 = 0.0
        r00=rerth/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
c
        blon = mod(blon,twopi)
        if(blon.lt.0. e 0) blon = twopi + blon
        r=r00*(tan(qtpi-is*blat/2))**cone
        x=    r*sin(cone*(blon-cenlon))
        y=-is*r*cos(cone*(blon-cenlon))
c
      endif
c
      cx = x
      cy = y
c
      return
      end
      subroutine makectls(igrd1,jgrd1,
     1                    rproj,rtruth,rorient,rdelx,rdely,
     2                    rcenlat,rcenlon,rlftgrd,rbtmgrd,nctl)
c
      dimension flat(igrd1*jgrd1),flon(igrd1*jgrd1),
     1          flatm(igrd1*jgrd1),iz(1)
c
      lngrd=igrd1*jgrd1
      call getgrd2(flat,flon,igrd1,jgrd1,
     1             rproj,rtruth,rorient,rdelx,rdely,
     2             rcenlat,rcenlon,rlftgrd,rbtmgrd)
c
      nctl=61
      open(nctl,file='polymask.ctl',status='unknown',form='formatted')
c
      ifh=0
      ihr=0
      iday=1
      imon=1
      iyr =99
      strlon=flon(1)*180./acos(-1.)
      endlon=flon(igrd1)*180./acos(-1.)
      if(endlon.lt.strlon) endlon=endlon+360.
      dlon=(endlon-strlon)/float(igrd1-1)
      print *,' strlon endlon dlon ',strlon,endlon,dlon
      do j=1,jgrd1
        i=(j-1)*igrd1+1
        flatm(j)=flat(i)*180./acos(-1.)
      enddo
      strlat=flatm(1)
      dlat=flatm(2)-flatm(1)
      iz(1)=0
      levr=1
      dx=rdelx/1000.
      call ctlhead(nctl,igrd1,jgrd1,levr,
     1             rproj,rlftgrd,rbtmgrd,rorient,dx,
     1             ihr,iday,imon,iyr,ifh,strlon,dlon,strlat,dlat,
     2             flatm,iz)
c
      write(nctl,100)
 100   format('vars 1')
      write(nctl,101)
 101  format('mask  0 99 mask')
      write(nctl,300)
 300  format('endvars')
c
      stop
      end
      subroutine ctlhead(nn,im,jm,km,proj,pi,pj,or,dx,
     1                  ihr,iday,imon,iyr,ifh,
     2                  strlon,dlon,strlat,dlat,
     3                  rlat,iz)
c$$$
      dimension rlat(jm),iz(km)
      character*2 hour,day
      character*3 mon(12)
      data mon/'jan','feb','mar','apr','may','jun',
     1         'jul','aug','sep','oct','nov','dec'/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      print *,' write control header to ctlhead '
      write(nn,101)
      write(nn,102)
      write(nn,103)
      write(nn,105)
      write(nn,106)
 101  format('dset ^polymask.data')
 102  format('options sequential')
 103  format('options byteswapped')
 105  format('undef -9.99e+33')
 106  format('title exp1')
c
      if( proj.eq.0.0 ) then
        write(nn,108) im,strlon,dlon
        write(nn,110) jm
        write(nn,111) (rlat(j),j=1,jm)
      else 
        write(nn,107) im,jm,pi,pj,or,dx
        imp=360.*111/dx
        strlonp=0.0
        dlonp=360./imp
        jmp=imp/4
        strlatp=0.0
        dlatp=dlonp
        write(nn,108) imp,strlonp,dlonp
        write(nn,109) jmp,strlatp,dlatp
      endif
 107  format('pdef',2i5,' nps',4f11.5)
 108  format('xdef',i5,' linear',2f9.3)
 109  format('ydef',i5,' linear',2f9.3)
 110  format('ydef',i5,' levels')
 111  format(10f8.3)
c
      write(nn,112) km
 112  format('zdef',i5,' levels 0 ')
      if( ihr.lt.10 ) then
        write(hour,90) ihr
      else
        write(hour,91) ihr
      endif
      if( iday.lt.10 ) then
        write(day,90) iday
      else
        write(day,91) iday
      endif
  90  format('0',i1)
  91  format(i2)
      if( ifh.eq.0 ) ifh=1
      write(nn,114) hour,day,mon(imon),iyr,ifh
 114  format('tdef 1 linear ',a2,'z',a2,a3,i2,i10,'hr')
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      subroutine getgrd2(flat,flon,igrd1,jgrd1,
     1                   rproj,rtruth,rorient,rdelx,rdely,
     2                   rcenlat,rcenlon,rlftgrd,rbtmgrd)
c
      dimension flat(igrd1*jgrd1),flon(igrd1*jgrd1)
c
      print 1234
 1234 format(' ==== in routine setgrd === ')
c
      hfpi = dasin(1.0d0)
      qtpi = hfpi * 0.5
      pi = 2.0 * hfpi
      twopi = 2.0 * pi
      rad = pi / 180.
      delx = rdelx
      dely = rdely
c
c --------- setup regional lat/lon and map factor -----
c
c if proj=0  do mercater projection
c if proj=1  do north polar projection
c if proj=-1 do south polar projection
c
      nproj = rproj
c
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
c ++++++++++++++++++++++++++++++++++++++
c polar projection
c ++++++++++++++++++++++++++++++++++++++
      truth  = rtruth * rad
      truth  = nproj * truth
      orient  = rorient * rad
      dlamda0 = orient + hfpi
      a2 =  6371200. * ( 1.0 + sin(truth) )
      radlat = rcenlat * rad
      radlon = rcenlon * rad - dlamda0
      radlat = nproj * radlat
      radlon = nproj * radlon
      yyy = a2 * cos( radlat )/(1. + sin( radlat ) )
      cenlat = rcenlat
      if( abs(cenlat) .eq. 90. ) then yyy = 0.0
      y00 = yyy * sin( radlon ) - ( rbtmgrd -1.) * dely
      x00 = yyy * cos( radlon ) - ( rlftgrd -1.) * delx
      print *,' delx x00 y00 ',delx,x00,y00
c
c =========
c           lat loop
      do 100 j = 1,jgrd1
      lats = j
      ijlats = (lats-1)*igrd1
      ys = y00 + (lats-1)*dely
c
      do 100 i=1,igrd1
      x = x00 + (i-1)*delx
      if( x .gt. 0. 0 ) then
         flons = atan(ys/x)
      else if ( x .lt. 0. 0 ) then
         flons = pi + atan(ys/x)
      else
         flons = hfpi
         if( ys .lt. 0. 0 ) flons = flons * 3.0
      endif
      flons = nproj * flons + dlamda0
      flons = mod(flons,twopi)
      if(flons.lt.0. 0) flons = twopi + flons
c
      rsoa2 = sqrt( x*x + ys*ys )/a2
      flats = hfpi - 2.0  * atan(rsoa2)
      flat(ijlats+i) = nproj * flats
      flon(ijlats+i) = flons
c
 100  continue
c
      else if ( nproj.eq.0 ) then
c
c ++++++++++++++++++++++++++++
c do mercater
c ++++++++++++++++++++++++++++
      truth  = rtruth * rad
      cenlat = rcenlat * rad
      cenlon = rcenlon * rad 
      a2 =  6371200. * cos( truth ) 
      x0 = 0.0
      y0 = a2 * log( abs( tan( qtpi + 0.5 * cenlat ) ) )
      x00 = - ( rlftgrd - 1.0 ) * delx
      y00 = - ( rbtmgrd - 1.0 ) * dely
      dlamda0 = 0.0
c
      do 200 j = 1,jgrd1
      lats = j
      ijlats = (lats-1)*igrd1
      ys = y00 + (lats-1)*dely + y0
c
       do 200 i=1,igrd1
         x = x00 + (i-1)*delx + x0
         flons = x / a2 + cenlon
c
         flats = 2.0 *( atan( exp( ys/a2 ) ) - qtpi )
         flat(ijlats+i) = flats
         flon(ijlats+i) = flons
c
 200  continue
      print *, 'flat ',(flat(i),i=1,10)
      print *, 'flon ',(flon(i),i=1,10)
c
      endif
c
c
      return
      end
