      program rsigdiff
c
      parameter (ntotal_=1)
      parameter(levmax=100)
c
      character*8 lab(4)
      character*8 labx(4)
c
      dimension idate(4),idatex(4)
c
      real, allocatable :: si(:),sl(:),
     1          flat(:), flon(:),
     1          fm2(:), fm2x(:), fm2y(:)
      real, allocatable :: q(:),te(:,:),
     1          uu(:,:),vv(:,:),rq(:,:)
      real, allocatable :: gz(:)
      real, allocatable :: dummy(:),ext(:)
c
      real, allocatable :: six(:),slx(:),
     1          flatx(:), flonx(:),
     1          fm2q(:), fm2qx(:), fm2qy(:)
      real, allocatable :: qx(:),tex(:,:),
     1          uux(:,:),vvx(:,:),rqx(:,:)
      real, allocatable :: gzx(:)
      real, allocatable :: dummyx(:),extx(:)
c
      character*128 file1,file2
	  namelist/namin/ file1,file2,levs,idim,jdim
      data n/11/
c
      read(5,namin)
c
      lngrd=idim*jdim
      km=levs
c
      allocate ( si(km+1),sl(km),
     1          flat(lngrd), flon(lngrd),
     1          fm2(lngrd), fm2x(lngrd), fm2y(lngrd) )
      allocate ( q(lngrd),te(lngrd,km),
     1          uu(lngrd,km),vv(lngrd,km),rq(lngrd,km*ntotal_))
      allocate ( gz(lngrd))
      allocate ( dummy(2*levmax),ext(512-(6+2*levmax)) )
c
      allocate ( six(km+1),slx(km),
     1          flatx(lngrd), flonx(lngrd),
     1          fm2q(lngrd), fm2qx(lngrd), fm2qy(lngrd) )
      allocate ( qx(lngrd),tex(lngrd,km),
     1          uux(lngrd,km),vvx(lngrd,km),rqx(lngrd,km*ntotal_) )
      allocate ( gzx(lngrd) )
      allocate ( dummyx(2*levmax),extx(512-(6+2*levmax)) )
c
c  first file
c
      open(unit=n,file=file1,status='unknown',
     1     form='unformatted')
c
      read(n) lab
      read(n) fhour,idate,(si(k),k=1,km+1),(sl(k),k=1,km)
     1         ,(dummy(i),i=1,2*levmax+1-km-1-km),ext
      read(n) (gz(i),i=1,lngrd)
      read(n) (q (i),i=1,lngrd)
      do k=1,km
      read(n) (te(i,k),i=1,lngrd)
      enddo
      do k=1,km
      read(n) (uu(i,k),i=1,lngrd)
      read(n) (vv(i,k),i=1,lngrd)
      enddo
      do k=1,km*ntotal_
      read(n) (rq(i,k),i=1,lngrd)
      enddo
      read(n) (fm2 (i),i=1,lngrd)
      read(n) (fm2x(i),i=1,lngrd)
      read(n) (fm2y(i),i=1,lngrd)
      read(n) (flat(i),i=1,lngrd)
      read(n) (flon(i),i=1,lngrd)
c
      open(unit=n,file=file2,status='unknown',
     1     form='unformatted')
      read(n) labx
      read(n) fhourx,idatex,(six(k),k=1,km+1),(slx(k),k=1,km)
     1         ,(dummyx(i),i=1,2*levmax+1-km-1-km),extx
      read(n) (gzx(i),i=1,lngrd)
      read(n) (qx (i),i=1,lngrd)
      do k=1,km
      read(n) (tex(i,k),i=1,lngrd)
      enddo
      do k=1,km
      read(n) (uux(i,k),i=1,lngrd)
      read(n) (vvx(i,k),i=1,lngrd)
      enddo
      do k=1,km*ntotal_
      read(n) (rqx(i,k),i=1,lngrd)
      enddo
      read(n) (fm2q (i),i=1,lngrd)
      read(n) (fm2qx(i),i=1,lngrd)
      read(n) (fm2qy(i),i=1,lngrd)
      read(n) (flatx(i),i=1,lngrd)
      read(n) (flonx(i),i=1,lngrd)
c
      print *,'lab  =',lab
      print *,'labx =',labx
      print *,'fhour =',fhour
      print *,'fhourx=',fhour
      print *,'idate =',idate
      print *,'idatex=',idate
      dsi=0.
      do k=1,km+1
       dsi=(si(k)-six(k))**2
      enddo
      print *,'dsi=',dsi
      dsl=0.
      do k=1,km
       dsl=(sl(k)-slx(k))**2
      enddo
      print *,'dsl=',dsl
      ddum=0.
      do i=1,2*levmax+1-km-1-km
       ddum=(dummy(i)-dummyx(i))**2
      enddo
      print *,'ddummy=',ddummy
      dext=0.
      do i=1,512-(6+2*levmax)
       dext=(ext(i)-extx(i))**2
      enddo
      print *,'dext=',dext
      dgz=0.
      dq =0.
      dfm2 =0.
      dfm2x=0.
      dfm2y=0.
      dflat=0.
      dflon=0.
      do i=1,lngrd
      dgz  =dgz  +(gz  (i)-  gzx(i))**2
      dq   =dq   +(q   (i)-   qx(i))**2
      dfm2 =dfm2 +(fm2 (i)-fm2q (i))**2
      dfm2x=dfm2x+(fm2x(i)-fm2qx(i))**2
      dfm2y=dfm2y+(fm2y(i)-fm2qy(i))**2
      dflat=dflat+(flat(i)-flatx(i))**2
      dflon=dflon+(flon(i)-flonx(i))**2
      enddo
      print *,'dgz,dq,dfm2,dfm2x,dfm2y,dflat,dflon=',
     1         dgz,dq,dfm2,dfm2x,dfm2y,dflat,dflon
      do k=1,km
       dte=0.
       do i=1,lngrd
       dte=dte+(te(i,k)-tex(i,k))**2
       enddo
       print *,'k,dte=',dte
      enddo
      do k=1,km
       dtu=0.
       dtv=0.
       do i=1,lngrd
       dtu=dtu+(uu(i,k)-uux(i,k))**2
       dtv=dtv+(vv(i,k)-vvx(i,k))**2
       enddo
       print *,'k,dtu,dtv=',dtu,dtv
      enddo
      do k=1,km*ntotal_
       dtr=0.
       do i=1,lngrd
       dtr=dtr+(rq(i,k)-rqx(i,k))**2
       enddo
       print *,'k,dtr=',dtr
      enddo
c
      stop
      end
