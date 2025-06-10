      program mapinfo
      implicit none
c#include <paramodel.h>
c#define DBG
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    mapino     Read related map information from grib header.
c                           Then calucalate the center position of the domain.
c                           All these will be saved in namelist &nameloc. The
c                           namelist is offered to rmtn to calculate the 
c                           orography in the specified regional domain.
c
c                           This program assumes that all the input standard
c                           pressure level fields are given in one grib format
c                           file (file name fn given as a namelist input).
c
c   prgmmr: masao kanamitsu date: 07-03-28
c           succeeded by kei yoshimura
c           haiqin li added regional to regional: 09-06-04
c
c attributes:
c   language: fortran 90.
c
      integer idim,jdim
c
      real rerth_
      parameter(rerth_=6.3712e+6)
      real g_
      parameter(g_=9.8000e+0)
      real pi_
      parameter(pi_=3.141593)

      integer nfld
      parameter(nfld=5)
c      integer mkd
c      parameter(mkd =(kdim*nfld+2)*2)
c
      real pi,twopi,hfpi,qtpi,rad
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.)
 
c
      real rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,
     1     rlftgrd,rbtmgrd,cproj,ctruth,corient,cdelx,cdely,
     2     ccenlat,ccenlon,clftgrd,cbtmgrd,rlat1,rlon1,rlat2,rlon2
      common /comrloc/
     1 rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd
     2,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
     3,rlat1,rlon1,rlat2,rlon2
      real  truth,orient,a2,x00,y1,y0,y00,cenlon,cenlat,rsoa2
      real  dlamda0
c
      real tensn
      parameter(tensn=10.0)
      real gamma,ps1,rh1,psk
      parameter(gamma=6.5e-3)
      integer lugb
      data lugb/10/
c
      integer iy,im,id,ih,indmois
      real    fh
      integer maxstdprs
      parameter(maxstdprs=100)
      real    stdprs(maxstdprs)
      character*128 fn
      data indmois/0/
      integer isfcp
      data isfcp/1/
      integer inst
      data inst/0/ !! 0: ncep rean grib id, 1: ecmwf rean grib id
      integer iintp
      data iintp/0/
      integer irhex
      data irhex/0/
      integer ihyb
      data ihyb/0/ !! 0: sigma level, 1: s-p hybrid level
c
      integer kdimprs
c
      real es,fpvs0,qs,xeps,xepsm1
c      parameter(xeps=rd_/rv_,xepsm1=rd_/rv_-1.)
c
      integer i,j,in,in2,k,n,m,lat,kdimi,kdimqi,ie,iw
c      real, allocatable :: out(:,:,:,:)
c
c
c   lugb : unit number for grib read file
c   fn   : input standard pressure level grid file name
c   iy,im,id,ih,fh   : input year,month,day,hour,forecast hour
c   stdprs : standard pressure level to be read in mb unit
c   kdimprs : number of standard pressure levels
c   indmois : =0 specific humidity ; =1 relative humidity
c   iorsmgsm : choice of output. =0 global spectral ; =1 regional
c   isfcp : =0 use p-level geoheight ; =1 use sfcP and topo
c   isfcp : =2 use sea level pressure and zero-topo
c   irhex : =0 use extrapolation for RH (default)
c   irhex : =1 do not extrapolate for upper than top level
c
      namelist/namp2s/fn,iy,im,id,ih,fh,stdprs,indmois,isfcp,inst,iintp
     $     ,irhex,ihyb
c
      do k=1,maxstdprs
        stdprs(k)=0.
      enddo
c
      read (5,namp2s)
      
      write(6,namp2s)
      
c
      if (ihyb.eq.1 .and. isfcp.ne.1) then
        print*,'if ihyb=1, isfcp should be 1'
        call abort
      endif
      kdimprs=1
      dowhile(stdprs(kdimprs).ne.0.)
        kdimprs=kdimprs+1
      enddo
      kdimprs=kdimprs-1
      print *,'kdimprs=',kdimprs
c
c      allocate (out(idimt,jdimhf,kdimprs,nfld))
c
      call prsgrbrd(lugb,fn,iy,im,id,ih,fh,stdprs,kdimprs,indmois,nfld,
     1              idim,jdim,isfcp,inst,iintp,ihyb)
      print*,'rproj,rorient,rtruth,rdelx,rdely=',
     1        rproj,rorient,rtruth,rdelx,rdely
c
c     to calculate rcenlat,rcenlat,rlftgrd,rbtmgrd (haiqin li)
c     now only support Mecator & Polar Stereographic projection
      if(rproj.eq.0) then ! Mecator projection  aaaaaa
        print*,'idim,jdim',idim,jdim
        rlftgrd=43
        rbtmgrd=35
        truth=rtruth*rad
        rlon1=rlon1*rad
        rlon1=mod(rlon1,twopi)
        rlat1=rlat1*rad
        a2=rerth_*cos(truth)
        x00=-(rlftgrd-1.)*rdelx - 0.5*rdelx
        cenlon=rlon1-x00/a2
        rcenlon=cenlon/rad
        if(rcenlon.gt.180) rcenlon=rcenlon-360.

        y1=(log(tan((rlat1/2.)+qtpi)))*a2
        y00=-(rbtmgrd-1.)*rdely-0.5*rdely
        y0=y1-y00
        cenlat=(atan(exp(y0/a2))-qtpi)*2.
        cenlat=mod(cenlat,pi_)
        rcenlat=cenlat/rad
        print*,'Mecator rcenlat,rcenlon',rcenlat,rcenlon
      elseif(abs(rproj).eq.1.) then ! Polar stereographic projection
        if(rproj.eq.1) then
           rcenlat=90.
           rtruth=60.
        elseif(rproj.eq.-1.) then
           rcenlat=-90.
           rtruth=-60.
        endif
        rcenlon=0.
        print*,'rcenlat,rcenlon,rorient',rcenlat,rcenlon,rorient
        rlat1=rlat1*rad
        rlon1=rlon1*rad
        truth  = rtruth * rad
        truth  = rproj * truth
        orient  = rorient * rad
        dlamda0 = orient + hfpi
        print*,'dlamda0=',dlamda0
        a2 =  rerth_ * ( 1.0 + sin(truth) )
        cenlat = rcenlat
        rlat1=rlat1/rproj
        rsoa2=tan((hfpi-rlat1)/2.0e0)
        rlon1=(rlon1-dlamda0)/rproj
        x00=((rsoa2*a2)**2/(1+(tan(rlon1-pi))**2))
        x00=-sqrt(x00)
!   (add a small number for integer truncation)
        y00=x00*tan(rlon1)-100.
        print*,'x00=',x00
        print*,'y00=',y00
        rlftgrd=1-x00/rdelx
        rbtmgrd=abs(1-y00/rdely)
        print*,'polar rlftgrd=',rlftgrd
        print*,'polar rbtmgrd=',rbtmgrd
      endif
c
        open(61,file='map.parm',form='formatted')
        write(61,*)'&NAMLOC'
        write(61,*)'RPROJ    = ',rproj,','
        write(61,*)'RTRUTH   = ',rtruth,','
        write(61,*)'RORIENT  = ',rorient,','   
        write(61,*)'RDELX    = ',rdelx,','
        write(61,*)'RDELY    = ',rdely,','
        write(61,*)'RCENLAT  = ',rcenlat,','
        write(61,*)'RCENLON  = ',rcenlon,','
        write(61,*)'RLFTGRD  = ',rlftgrd,','
        write(61,*)'RBTMGRD  = ',rbtmgrd,','
        write(61,*)'CPROJ    = 0.0,'
        write(61,*)'CTRUTH   = 0.0,'
        write(61,*)'CORIENT  = 0.0,'
        write(61,*)'CDELX    = 0.0,'
        write(61,*)'CDELY    = 0.0,'
        write(61,*)'CCENLAT  = 0.0,'
        write(61,*)'CCENLON  = 0.0,'
        write(61,*)'CLFTGRD  = 0.0,'
        write(61,*)'CBTMGRD  = 0.0,'
        write(61,*)'/'
c

      end

      subroutine prsgrbrd(lugb,fn,iy,im,id,ih,fh,
     1                    stdprs,kdimprs,indmois,nfld,
     2                    idim,jdim,isfcp,inst,iintp,ihyb)
      implicit none
      real rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,
     1     rlftgrd,rbtmgrd,cproj,ctruth,corient,cdelx,cdely,
     2     ccenlat,ccenlon,clftgrd,cbtmgrd,rlat1,rlon1,rlat2,rlon2
      common /comrloc/
     1 rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd
     2,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
     3,rlat1,rlon1,rlat2,rlon2
c
c  read standard pressure level grib file, interpolates to model grid.
c
c   abs(proj)=0  : mercator
c            =1  : polar stereopgraphic
c            =2  : lambert
c            =3  : lat/lon
c            =4  : gaussian
c    negative value indicates origin in southern hemisphere
c
      real pi_
      parameter(pi_=3.141593)
c
      integer       lugb
      character*128 fn
      integer       idim,jdim
      integer       iy,im,id,ih
      real          fh
      integer       kdimprs
      real          stdprs(kdimprs)
      integer       indmois
      integer       nfld
      integer       isfcp
      integer       inst
      integer       iintp
      integer       ihyb
c
      real, allocatable :: data(:)
c
      integer,allocatable ::  kpd(:,:,:)
      integer k,kgau,i,ij,kpds5,kgds1,ijmdim,imaxgrb,jmaxgrb,ijdim
      real cotru
      real dlon,dlat,wlon,rnlat
      integer iret
      real,allocatable :: dummy1(:,:)
      integer ij2,ii,jj,j
c
c  max possible number of grid points on input grib file
c
      integer imdata,jmdata
      parameter(imdata=720,jmdata=361)
c
      integer n
      integer kpds5u,kpds5v,kpds5t,kpds5z,kpds5q,kpds5rh,kpds5p
      integer kpds5slp
#ifdef DBG
      print *,'lugb,fn=',lugb,fn
#endif
c
c
      if(fn(1:4).eq.'    ') then
         print *,'pressure grib file name empty.'
         call abort
      endif
c
      allocate (kpd(25,kdimprs,nfld))
c
      ijmdim=imdata*jmdata
      print*,'1'
      allocate (data(ijmdim))
c
c  set kpd array
c ncep rean
c  n=1; u  (kpds5=33)
c  n=2: v  (kpds5=34)
c  n=3: t  (kpds5=11)
c  n=4: q  (kpds5=51) when indmois=0
c  n=4: rh (kpds5=52) when indmois=1
c  n=5: z  (kpds5= 7) when isfcp=0
c  n=5: ps (kpds5= 1) when isfcp=1
c  n=5: slp(kpds5= 2) when isfcp=2 (use zero-topo)
c ecmwf rean
c  n=1; u  (kpds5=131)
c  n=2: v  (kpds5=132)
c  n=3: t  (kpds5=130)
c  n=4: q  (kpds5=133) when indmois=0
c  n=4: rh (kpds5=157) when indmois=1
c  n=5: z  (kpds5=129) when isfcp=0 ; be aware of unit [m2/s2]
c  n=5: ps (kpds5= 1) when isfcp=1
c
c   kpds6=100 for isobaric level
c   kpds7=isobaric level pressure in hPa
c
      if (inst.eq.0) then !! ncep rean
         kpds5u=33
         kpds5v=34
         kpds5t=11
         kpds5q=51
         kpds5rh=52
         kpds5z=7
         kpds5p=1
         kpds5slp=2
      elseif (inst.eq.1) then !! ecmwf rean
         kpds5u=131
         kpds5v=132
         kpds5t=130
         kpds5q=133
         kpds5rh=157
         kpds5z=129
         kpds5p=1
         kpds5slp=2
      endif

      do n=1,nfld
        do k=1,kdimprs
          do i=1,25
            kpd(i,k,n)=-1
          enddo
        enddo
      enddo
      do k=1,kdimprs
        do n=1,nfld
          if (ihyb.eq.0) then
            kpd(6,k,n)=100
          else
            kpd(6,k,n)=109
          endif
          kpd(7,k,n)=stdprs(k)
        enddo
        kpd(5,k,1)=kpds5u
        kpd(5,k,2)=kpds5v
        kpd(5,k,3)=kpds5t
        if(indmois.eq.0) then
          kpd(5,k,4)=kpds5q
        elseif(indmois.eq.1) then
          kpd(5,k,4)=kpds5rh
        else
          print *,'wrong indmois'
          call abort
        endif
        if (isfcp.eq.0) then
          kpd(5,k,5)= kpds5z
        elseif (isfcp.eq.1) then
          if (k.eq.1) then !! surface pressure
            kpd(5,k,5)=kpds5p
            kpd(6,k,5)=1
            kpd(7,k,5)=0
          elseif (k.eq.2) then !! topography
            kpd(5,k,5)=kpds5z
            kpd(6,k,5)=1
            kpd(7,k,5)=0
          else
            kpd(5,k,5)=kpds5t      !! temporary
          endif
        elseif (isfcp.eq.2) then
          if (k.eq.1) then !! sealevel pressure
            kpd(5,k,5)=kpds5slp
            kpd(6,k,5)=102
            kpd(7,k,5)=0
          else
            kpd(5,k,5)=kpds5t      !! temporary
          endif
        endif
      enddo
c      kpd(5,2,5)=1
c      kpd(6,2,5)=1
c      kpd(7,2,5)=0
c      kpd(5,3,5)=7
c      kpd(6,3,5)=1
c      kpd(7,3,5)=0
c
c  read grib file
c
      do n=1,nfld
c         print*,'in read grib n=',n
        do k=1,kdimprs
 82       continue
          print *,'== start reading ',kpd(5,k,n),kpd(6,k,n),kpd(7,k,n)
          call fixrd2(lugb,fn,kpd(1,1,1),    !  kpd(1,k,n),
     &               iy,im,id,ih,fh,
     &               imdata,jmdata,
     &               imaxgrb,jmaxgrb,
     &               kpds5,kgds1,
     1               idim,jdim,
     &               rproj,rorient,rtruth,cotru,
     &               rdelx,rdely,rlat1,rlat2,rlon1,rlon2,
     &               data,iret)
          if(iret.gt.0) then
            if (isfcp.eq.2.and.kpd(5,k,n).eq.2) then
                                !! allow either kpd5=1 or 2 for mslp
              kpd(5,k,n)=1
              print *,'== instead, look for '
     $             ,kpd(5,k,n),kpd(6,k,n),kpd(7,k,n)
              goto 82
            else
              print *,kpd(5,k,n),kpd(6,k,n),kpd(7,k,n),' read failed'
              call abort
            endif
          endif
        enddo
      enddo
c
      deallocate (data)
      deallocate (kpd)

      return
      endsubroutine prsgrbrd

      subroutine fixrd2(lugb,fngrib,kkpds,
     &                 iy,im,id,ih,fh,imdata,jmdata,
     &                 imax,jmax,kkpds5,kkgds1,
     1                 idim,jdim,
     &                 proji,orienti,truthi,cotrui,
     &                 delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i,
#ifdef G2G
     &                 dlon,dlat,wlon,rnlat,
#endif
     &                 data,iret)
c#include <paramodel.h>
#define REAL4_W3LIB
#define DBG
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c read in grib climatology files.
c
c interpolate climatology to the dates
c
c grib file should allow all the necessary parameters to be extracted fr
c the description records.
c
c
c  nvalid:  analysis later than (current date - nvalid) is regarded as
c           valid for current analysis
c
#ifdef ADJSOILW
      parameter(nvalid=1)
#else
      parameter(nvalid=5)
#endif
c
c  iret=0  ... successfully read
c      =-1 ... record with unmatched dates exist
c
      character*128 fngrib
      integer iy,im,id,ih
c
      real data(imdata*jmdata)
      real, allocatable :: data_next(:)
c
      logical*1 lbms(imdata*jmdata)
c
      parameter(mbuf=1024*128*64)
      character*1 cbuf(mbuf)
c
      integer kpds(25),kgds(22),kens(5)
      integer jpds(25),jgds(22),jens(5)
      integer kkpds(25),kpds0(25)
c
      integer idim,jdim
c
      character*80 asgnstr
c
      data msk1/3840000/,msk2/2400000/
c
#ifdef REAL4_W3LIB
      integer*4 lugb4,msk14,msk24,mnum4,mbuf4
      integer*4 nlen4,nnum4,iret4
      integer*4 ndata4
      real*4, allocatable :: data4(:)
      integer*4 lskip4,lgrib4,lret4
      integer*4 n4,jpds4(25),jgds4(22),jens4(5)
      integer*4 k4,kpds4(25),kgds4(22),kens4(5)
#endif
c
#ifdef FASTBAREAD
      character, allocatable :: bbuf(:)
#endif
c
      logical lclim
c
c julian day of the middle of each month
c
      real dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     1           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
      save dayhf
c
#ifdef REAL4_W3LIB
      allocate (data4(imdata*jmdata))
#endif
c
      print*,'bf iy,im,id,ih,fh',iy,im,id,ih,fh
      call vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
      print*,'af iy,im,id,ih,fh',iy,im,id,ih,fh
c
      close(lugb)
      call numchar(fngrib,nch)
#ifdef ASSIGN
      if(lugb.lt.10) then
        write(asgnstr,'(22hassign -s unblocked u:,i1)') lugb
      else
        write(asgnstr,'(22hassign -s unblocked u:,i2)') lugb
      endif
      call assign('assign -R')
      call assign(asgnstr)
      open(unit=lugb,file=fngrib(1:nch),status='old',form='unformatted',
     1     err=910)
#else
#ifdef FASTBAREAD
      call nainit(lugb,fngrib(1:nch),isize,iret)
      allocate (bbuf(isize))
      call naopen(lugb,fngrib(1:nch),bbuf,isize,iret)
#else
      call baopen(lugb,fngrib(1:nch),iret)
#endif
      if(iret.ne.0) go to 910
#endif
      go to 911
  910 continue
      print *,'error in opening file ',fngrib(1:nch)
      call abort
c
  911 continue
      write(6,*) ' file ',fngrib(1:nch),' opened. unit=',lugb
c
c  get grib index buffer
c
      mnum=0
#ifdef REAL4_W3LIB
      lugb4=lugb
      msk14=msk1
      msk24=msk2
      mnum4=mnum
      mbuf4=mbuf
#ifdef FASTBAREAD
      call ngetgir(bbuf,isize,msk14,msk24,mnum4,mbuf4,
     1            cbuf,nlen4,nnum4,iret4)
#else
      call getgir(lugb4,msk14,msk24,mnum4,mbuf4,
     1            cbuf,nlen4,nnum4,iret4)
#endif
      nlen=nlen4
      nnum=nnum4
      iret=iret4
#else
#ifdef FASTBAREAD
      call ngetgir(bbuf,isize,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#else
      call getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#endif
#endif
#ifdef DBG
      print *,'nlen=',nlen,' nnum=',nnum
#endif
      if(iret.ne.0) then
        print *,'error.  cbuf length too short in ngetgir'
        call abort
      endif
      if(nnum.eq.0) then
        print *,'error.  not a grib file. detected in ngetgir'
        call abort
      endif
      if(nlen.eq.0) then
        print *,'error.  nlen=0.  detected in ngetgir'
        call abort
      endif
      if(nlen.gt.imdata*jmdata) then
        print *,'error.  nlen .gt. imdata*jmdata  detected in ngetgir'
        call abort
      endif
c
c  find file type climatology or analysis
c
      do i=1,25
        jpds(i)=-1
      enddo
      do i=1,22
        jgds(i)=-1
      enddo
      do i=1,5
        jens(i)=-1
      enddo
      do i=1,25
        jpds(i)=kkpds(i)
      enddo
c
c fix for ecmwf grib
c
      if ( (jpds(5).eq.11.or.jpds(5).eq.91) .and.
     &     (jpds(7).eq.-1.or.jpds(7).eq.0) ) then
        jpds(5)=-1
        jpds(7)=0
        write(6,*) 'Changing jpds5 and jpds7 to allow ',
     1             'for ECMWF and NCEP ice and sst'
      endif
c
#ifdef DBG
      write(6,*) ' Searching the following fields with kpds=.'
      write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
      write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
      write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
      n=0
#ifdef REAL4_W3LIB
      nlen4=nlen
      nnum4=nnum
      n4=n
      do i=1,25
        jpds4(i)=jpds(i)
      enddo
      do i=1,22
        jgds4(i)=jgds(i)
      enddo
      do i=1,5
        jens4(i)=jens(i)
      enddo
      call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &             k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
      k=k4
      do i=1,25
        kpds(i)=kpds4(i)
      enddo
      do i=1,22
        kgds(i)=kgds4(i)
      enddo
      do i=1,5
        kens(i)=kens4(i)
      enddo
      lskip=lskip4
      lgrib=lgrib4
      iret=iret4
#else
      call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &             k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
#ifdef DBG
      write(6,*) 'first grib record in the  grib file',fngrib(1:nch)
      write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
      write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
      write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
c
      if(lgrib.eq.0) then
        write(6,*) ' error in getgbss.  No matching records.'
        iret=99
        return
c        call abort
      endif
c
      do i=1,25
        kpds0(i)=kpds(i)
      enddo
      kpds0(4)=-1
      kpds0(18)=-1
c
clu_rev: manipulate kpds(13) and kpds(15) for fixed fields
c
c neet to clean up  mk
c
      if(kpds(5).eq.236 .or. kpds(5).eq.255) then
         kpds(13) = 4
         kpds(15) = 1
         kpds(16) = 51
         print *,' quick fix for noah fixed fields'
      endif
c
      if(kpds(16).eq.51) then
        write(6,*) ' climatology file.'
        lclim=.true.
      else
        write(6,*) ' analysis file.'
        lclim=.false.
      endif
c
c   handling climatology file
c
      if(lclim) then
c
c       find average type
c         weekly,biweekly,monthly,seasonal,annual
c
c  kpds(13)=4 & kpds(15)=1 .. annual mean
c  kpds(13)=3 & kpds(15)=1 .. monthly mean
c  kpds(13)=2 & kpds(15)=7 .. weekly mean
c  kpds(13)=2 & kpds(15)=14.. bi-weekly mean
c
        if(kpds(13).eq.2.and.kpds(15).eq.7) then
          write(6,*) ' this is weekly mean climatology'
          write(6,*) ' cannot process.'
          call abort
        elseif(kpds(13).eq.2.and.kpds(15).eq.14) then
          write(6,*) ' this is bi-weekly mean climatology'
          write(6,*) ' cannot process.'
          call abort
        elseif(kpds(13).eq.3.and.kpds(15).le.1) then
          write(6,*) ' this is monthly mean climatology'
          monend=12
          do mm=1,monend
            mmm=mm
            mmp=mm+1
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
              mon1=mmm
              mon2=mmp
              go to 20
            endif
          enddo
          print *,'wrong rjday',rjday
          call abort
 20       continue
          ijmax=0
          do nn=1,2
            if(nn.eq.2) then
              allocate (data_next(ijmax))
            endif
            lskip = -1
            n=0
            do i=1,25
              jpds(i)=kpds0(i)
            enddo
            jpds(24)=-1
            jpds(25)=-1
            if(nn.eq.1) jpds( 9)=mon1
            if(nn.eq.2) jpds( 9)=mon2
            if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
            write(6,*) ' Searching the following fields with kpds=.'
            write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
            write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
            write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
            nlen4=nlen
            nnum4=nnum
            n4=n
            do i=1,25
              jpds4(i)=jpds(i)
            enddo
            do i=1,22
              jgds4(i)=jgds(i)
            enddo
            do i=1,5
              jens4(i)=jens(i)
            enddo
            call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &                   k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
            k=k4
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            do i=1,5
              kens(i)=kens4(i)
            enddo
            lskip=lskip4
            lgrib=lgrib4
            iret=iret4
#else
            call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &                   k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
            if(lgrib.eq.0) then
#ifdef DBG
              write(6,*) ' No matching record found.'
              write(6,*) ' The last header records read:'
              write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
              write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
              write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
              call abort
            endif
c
            write(6,*) ' Matching record found for mon=',jpds(9)
#ifdef DBG
            write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
            write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
            write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
#endif
            do i=1,imdata*jmdata
              lbms(i)=.true.
            enddo
#ifdef REAL4_W3LIB
            lugb4=lugb
#ifdef FASTBAREAD
            call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#else
            call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#endif
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            ndata=ndata4
            if(nn.eq.1) then
              do i=1,ndata
                data(i)=data4(i)
              enddo
            else
              do i=1,ndata
                data_next(i)=data4(i)
              enddo
            endif
#else
            if(nn.eq.1) then
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#endif
            else
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#endif
            endif
#endif
            if(ndata.eq.0) then
              write(6,*) ' error in nrdgb'
              call abort
            endif
            imax=kgds(2)
            jmax=kgds(3)
            ijmax=imax*jmax
#ifdef DBG
            write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
          enddo
c
          wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
#ifdef DBG
          print *,'rjday,mon1,mon2,wei1,wei2=',
     1             rjday,mon1,mon2,wei1,wei2
#endif
          do i=1,ijmax
            data(i)=wei1*data(i)+wei2*data_next(i)
          enddo
          deallocate (data_next)
        elseif(kpds(13).eq.4.and.kpds(15).eq.3) then
          write(6,*) ' this is seasonal mean climatology'
          monend=4
          is=im/3+1
          if(is.eq.5) is=1
          is1=mon1/3+1
          if(is1.eq.5) is1=1
          is2=mon2/3+1
          if(is2.eq.5) is2=1
          do mm=1,monend
            mmm=mm*3-2
            mmp=(mm+1)*3-2
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
              mon1=mmm
              mon2=mmp
              go to 30
            endif
          enddo
          print *,'wrong rjday',rjday
          call abort
   30     continue
          ijmax=0
          do nn=1,2
            if(nn.eq.2) then
              allocate (data_next(ijmax))
            endif
            do i=1,25
              jpds(i)=kpds0(i)
            enddo
            jpds(24)=-1
            jpds(25)=-1
            n=0
            if(nn.eq.1) then
              isx=is1
            else
              isx=is2
            endif
            if(isx.eq.1) jpds(9)=12
            if(isx.eq.2) jpds(9)=3
            if(isx.eq.3) jpds(9)=6
            if(isx.eq.4) jpds(9)=9
            if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
            write(6,*) ' Searching the following fields with kpds=.'
            write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
            write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
            write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
            nlen4=nlen
            nnum4=nnum
            n4=n
            do i=1,25
              jpds4(i)=jpds(i)
            enddo
            do i=1,22
              jgds4(i)=jgds(i)
            enddo
            do i=1,5
              jens4(i)=jens(i)
            enddo
            call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &                   k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
            k=k4
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            do i=1,5
              kens(i)=kens4(i)
            enddo
            lskip=lskip4
            lgrib=lgrib4
            iret=iret4
#else
            call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &                   k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
            if(lgrib.eq.0) then
#ifdef DBG
              write(6,*) ' No matching record found.'
              write(6,*) ' The last header records read:'
              write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
              write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
              write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
              call abort
            endif
            do i=1,imdata*jmdata
              lbms(i)=.true.
            enddo
#ifdef REAL4_W3LIB
            lugb4=lugb
#ifdef FASTBAREAD
            call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#else
            call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#endif
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            ndata=ndata4
            if(nn.eq.1) then
              do i=1,ndata
                data(i)=data4(i)
              enddo
            else
              do i=1,ndata
                data_next(i)=data4(i)
              enddo
            endif
#else
            if(nn.eq.1) then
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#endif
            else
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#endif
            endif
#endif
            if(ndata.eq.0) then
              write(6,*) ' error in nrdgb'
              call abort
            endif
            imax=kgds(2)
            jmax=kgds(3)
            ijmax=imax*jmax
          enddo
          wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
#ifdef DBG
          print *,'rjday=',rjday
          print *,'mon1 =',mon1 ,' mon2=',mon2
          print *,'wei1 =',wei1 ,' wei2=',wei2
          print *,'ses1 =', is1 ,' ses2=', is2
#endif
          do i=1,ijmax
            data(i)=wei1*data(i)+wei2*data_next(i)
          enddo
          deallocate (data_next)
        elseif(kpds(13).eq.4.and.kpds(15).eq.1) then
          write(6,*) ' this is annual mean climatology'
          monend=-1
#ifdef REAL4_W3LIB
          lugb4=lugb
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2              data4,6)
#else
          call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2              data4,6)
#endif
          do i=1,25
            kpds(i)=kpds4(i)
          enddo
          do i=1,22
            kgds(i)=kgds4(i)
          enddo
          ndata=ndata4
          do i=1,ndata
            data(i)=data4(i)
          enddo
#else
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
          call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
          if(ndata.eq.0) then
            write(6,*) ' error in nrdgb'
            write(6,*) ' kpds=',kpds
            write(6,*) ' kgds=',kgds
            call abort
          endif
          imax=kgds(2)
          jmax=kgds(3)
          ijmax=imax*jmax
        else
          write(6,*) ' climatology file average period unknown.'
          write(6,*) ' kpds(13)=',kpds(13),' kpds(15)=',kpds(15)
          call abort
        endif
      else
c-------------------------------------------------------------------
c
c  handling analysis file
c
c  find record for the given hour/day/month/year
c
        monend=0
        nrept=0
        do i=1,25
          kpds(i)=kpds0(i)
        enddo
        iyr=jy
        imo=jm
!        idy=jd
!       debug haiqin li
        idy=id
        ihr=jh
!   50   continue
        jpds( 8)=iyr
        jpds( 9)=imo
        jpds(10)=idy
        if(ihr.eq.-1) then
c
c t.r.    jpds(11)=0
          jpds(11)=-1
        else
c
c          jpds(11)=ihr
!        debug haiqin li
          jpds(11)=0
c          jpds(11)=-1
        endif
#ifdef DBG
        write(6,*) ' analysis grib file'
        write(6,*) ' Searching the following fields with kpds=.'
        write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
        write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
        write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
        n=0
#ifdef REAL4_W3LIB
        nlen4=nlen
        nnum4=nnum
        n4=n
        do i=1,25
          jpds4(i)=jpds(i)
        enddo
        do i=1,22
          jgds4(i)=jgds(i)
        enddo
        do i=1,5
          jens4(i)=jens(i)
        enddo
        call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &               k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
        k=k4
        do i=1,25
          kpds(i)=kpds4(i)
        enddo
        do i=1,22
          kgds(i)=kgds4(i)
        enddo
        do i=1,5
          kens(i)=kens4(i)
        enddo
        lskip=lskip4
        lgrib=lgrib4
        iret=iret4
#else
        call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &               k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
        if(lgrib.ne.0) then
          if(nrept.le.iabs(nvalid).and.nrept.ne.0) then
            write(6,*) '<warning:cycl> grib record date does not match',
     1                 ' but within ',iabs(nvalid),' days.'
            write(6,*) '<warning:cycl> assume it is valid.'
          endif
          write(6,*) ' Matching record found'
#ifdef DBG
          write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
          write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
          write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
          write(6,*) ' Now reading grib record'
#endif
#ifdef REAL4_W3LIB
          lugb4=lugb
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     1              data4,6)
#else
          call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,data4,6)
#endif
          do i=1,25
            kpds(i)=kpds4(i)
          enddo
          do i=1,22
            kgds(i)=kgds4(i)
          enddo
          ndata=ndata4
          do i=1,ndata
            data(i)=data4(i)
          enddo
#else
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
          call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
          if(ndata.eq.0) then
            write(6,*) ' error in nrdgb'
            write(6,*) ' kpds=',kpds
            write(6,*) ' kgds=',kgds
            write(6,*) ' lgrib,lskip=',lgrib,lskip
            call abort
          endif
#ifdef DBG
          write(6,*) ' matching grib record successfully read.'
#endif
          imax=kgds(2)
          jmax=kgds(3)
          ijmax=imax*jmax
        else
           write(6,*) ' no matching dates found. abort.'
            call abort
c          endif
        endif
      endif
c
c  finished reading grib file.
c
c  substitution to the output array
c
   80 continue
      call rmaxmin(data,imax*jmax,1,1,1,'input data')
c
#ifdef DBG
c     if(kkpds(5).eq.225) then
c       print *,'input uninterpolated vegtype'
c       call nntprt(data,imax,jmax)
c     endif
c     if(kkpds(5).eq.87) then
c       print *,'input uninterpolated vegfrac'
c       call nntprt(data,imax,jmax)
c     endif
c     if(kkpds(5).eq.230) then
c       print *,'input uninterpolated soiltype'
c       call nntprt(data,imax,jmax)
c     endif
#endif
c
#ifdef DBG
      write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
      proji=0.
      orienti=0.
      truthi=0.
      cotrui=0.
      call rgetarea(kgds,proji,orienti,truthi,cotrui,
     1              delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i)
#ifdef DBG
      write(6,*) 'proji,orienti,truthi,cotrui='
      write(6,*)  proji,orienti,truthi,cotrui
      write(6,*) 'delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i='
      write(6,*)  delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i
#endif
!      call getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr)
!#ifdef DBG
!      write(6,*) 'dlat,dlon,rslat,rnlat,wlon,elon,ijordr ='
!      write(6,*)  dlat,dlon,rslat,rnlat,wlon,elon,ijordr
!#endif
!      call subst(data,imax,jmax,dlon,dlat,ijordr)
c
      write(6,*) 'imax,jmax of grib=',imax,jmax
c
      kkpds5=kkpds(5)
      kkgds1=kgds(1)
      idim=kgds(2)
      jdim=kgds(3)
c
#ifdef REAL4_W3LIB
      deallocate (data4)
#endif
#ifdef FASTBAREAD
      deallocate (bbuf)
#endif
      close(lugb)
c
      return
      endsubroutine fixrd2

      subroutine vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
c#include <paramodel.h>
c
c julian day of the middle of each month
c
      dimension dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     1           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
      save dayhf
c
c number of days in a month
c
      dimension mjday(12)
      data mjday/31,28,31,30,31,30,31,31,30,31,30,31/
      save mjday
c
c julian day of the first day of the month
c
      dimension fjday(12)
c
      logical ijordr
c
      iret=0
c
      monend=9999
c
      fjday(1)=1.
      do mon=2,12
        fjday(mon)=fjday(mon-1)+float(mjday(mon-1))
      enddo
c
c  get julian day of the iy/im/id/ih provided.
c
      rjday=0.
      imm=im-1
      if(imm.ge.1) then
        do mon=1,imm
          rjday=rjday+mjday(mon)
        enddo
      endif
      rjday=rjday+id+float(ih)/24.0
      rjday=rjday+fh/24.0
      rjday = mod(rjday,365.)
      if(rjday.eq.0.) rjday = 365.
      if(rjday.le.0..or.rjday.gt.365.) then
        print *,'wrong rjday',rjday
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
      do mon=1,11
        if(rjday.lt.fjday(mon+1)) then
          mfmon=mon
          go to 10
        endif
      enddo
      mfmon=12
   10 continue
      if(rjday.lt.dayhf(1)) rjday=rjday+365.
c
c  compute jy,jm,jd,jh of forecast
c
      jy=iy
      jm=im
      jd=id
      incdy=int(fh/24.)
      jh=ih+mod(fh,24.)
      incdy=incdy+jh/24
      jh=mod(jh,24)
      if(incdy.ge.1) then
      do incd=1,incdy
        jd=jd+1
        if(jm.eq.4.or.jm.eq.6.or.jm.eq.9.or.jm.eq.11) then
          mondy=30
        elseif(jm.eq.2) then
          if(mod(jy,4).eq.0) then
            mondy=29
          else
            mondy=28
          endif
        else
          mondy=31
        endif
        if(jd.gt.mondy) then
          jm=jm+1
          jd=1
          if(jm.gt.12) then
            jy=jy+1
            jm=1
          endif
        endif
      enddo
      endif
      write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
c
      return
      endsubroutine vday

      subroutine rgetarea(kgds,proj,orient,truth,cotru,
     1                   delx,dely,rlat1,rlat2,rlon1,rlon2)
c#include <paramodel.h>
c
c  get area of the grib record
c
cfpp$ noconcur r
      dimension kgds(25)
      logical ijordr
c
      real rerth_
      parameter(rerth_=6.3712e+6)
      real pi_
      parameter(pi_=3.141593)
      rad=pi_/180.
#ifdef DBG
      write(6,*) ' kgds( 1-10)=',(kgds(j),j= 1,10)
      write(6,*) ' kgds(11-20)=',(kgds(j),j=11,20)
      write(6,*) ' kgds(21-  )=',(kgds(j),j=21,25)
#endif
c
      idrt=kgds(1)
      rlat1=kgds(4)*1.e-3
      rlon1=kgds(5)*1.e-3

          if(idrt.eq.0) then
            proj=3
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point
            delx=kgds(09)  ! dx (milledegree) on truth latitude
            dely=kgds(10)  ! dy (milledegree) on truth latitude
c           orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
            if(kgds(11).eq.0.) proj = -3.
          endif
          if( idrt.eq.4 ) then      ! gaussian projection
            proj=4
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point
            delx=kgds(09)  ! dx (milledegree) on truth latitude
            dely=kgds(10)  ! dy (milledegree) on truth latitude
c           orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
            if(kgds(11).eq.0.) proj = -4.
          endif
          if( idrt.eq.1 ) then      ! mercater projection
            proj=0
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point

c regional mercator projection that crosses 0E.
c hyodae/kana 6/3/2008
              if(rlon1.gt.0.and.rlon1>rlon2) then
                  rlon2=rlon2+360.
               endif

c            delx=(rlon2-rlon1)/float(kgds(02)-1)
            delx=kgds(12)  ! dx (meter) on truth latitude
            dely=kgds(13)  ! dy (meter) on truth latitude
            orient=rlon1
            truth=kgds(09)*1.e-3 ! truth latitude
            cotru=truth          ! co-truth latitude
c            delx=delx*cos(truth*rad)*rerth_*rad
          endif
          if( idrt.eq.5 ) then      ! polar projection
            truth=60.0              ! truth latitude
            cotru=60.0              ! co-truth latitude
            orient=kgds(07)*1.e-3   ! orientation
            delx=kgds(08)  ! dx (meter) on 60 deg
            dely=kgds(09)  ! dy (meter) on 60 deg
            iproj=kgds(10) ! polar projection (first bit 0:north;1=south)
            if( iproj.eq.0 ) then
              proj=1.0
            else
              proj=-1.0
            endif
          endif
          if( idrt.eq.3 ) then      ! lambert projection
            orient=kgds(07)*1.e-3   ! orientation
            delx=kgds(08)  ! dx (meter) on 60 deg
            dely=kgds(09)  ! dy (meter) on 60 deg
            iproj=kgds(10) !  projection (first bit 0:north;1=south)
            if( iproj.eq.0 ) then
              proj=2.0
            else
              proj=-2.0
            endif
            truth=kgds(12)*1.e-3 ! the 1st lat from pole to cut
            cotru=kgds(13)*1.e-3 ! the 2nd lat from pole to cut
          endif
c
      return
      endsubroutine rgetarea

      subroutine numchar(char,num)
c
      implicit none
      character*128 char
      integer num,n
c
      num=1
      dowhile(char(num:num).ne.' ')
        num=num+1
      enddo
      num=num-1
      do n=num+1,128
        char(n:n)=' '
      enddo
      return
      endsubroutine numchar



