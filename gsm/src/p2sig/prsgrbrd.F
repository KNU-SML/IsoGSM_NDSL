#define PI pi
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
      subroutine prsgrbrd(lugb,fn,iy,im,id,ih,fh,
     1                    stdprs,kdimprs,indmois,nfld,
     2                    out,idim,jdim,isfcp,inst,iintp,ihyb,izrtr)
      implicit none
#include <paramodel.h>
#define G2G
c
c  read standard pressure level grib file, interpolates to model grid.
c
c  input:
c
c    lugb    int                   .. unit number for reading grib
c    fn      char                  .. grib file name
c    idim    int                   .. output x-dimension (model x-dimension)
c    jdim    int                   .. output y-dimension (model y-dimension)
c    iy      int                   .. 4-digit year
c    im      int                   .. month
c    id      int                   .. day
c    ih      int                   .. hour
c    fh      real                  .. forecast hour
c    stdprs  real                  .. standard pressure levels to extract
c    kdimprs int                   .. number of standard pressure levels
c    indmois int                   .. moisture variable indicator :: 0=specific humidity
c    nfld    int                   .. number of dependent variables (normally 5)
c                                                                    1=relative humidity
c
c  output:
c   
c    out    real array (idim*jdim,kdimprs) .. output field
c
c  the map related parameters appear in this program are for grib
c  files and not those of the forecast program
c
c  note on map projection:
c
c   abs(proj)=0  : mercator
c            =1  : polar stereopgraphic
c            =2  : lambert
c            =3  : lat/lon
c            =4  : gaussian
c    negative value indicates origin in southern hemisphere
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
      real          out(idim*jdim,kdimprs,*)
      integer       isfcp
      integer       inst
      integer       iintp
      integer       ihyb
      integer       izrtr
c
      real, allocatable :: data(:)
c
      integer,allocatable ::  kpd(:,:,:)
      integer k,kgau,i,ij,kpds5,kgds1,ijmdim,imaxgrb,jmaxgrb,ijdim
#ifdef G2G
      real proj,orient,truth,cotru
      real delx,dely
      real rlat1,rlat2,rlon1,rlon2
#else
#include <rscomloc.h>
      real cotru
#endif
      real dlon,dlat,wlon,rnlat
      integer iret
      real,allocatable :: dummy1(:,:)
      real dummy2(idim,jdim)
      integer ij2,ii,jj,j
      real          tmp(idim*jdim)
c
c  max possible number of grid points on input grib file
c
      integer imdata,jmdata
      parameter(imdata=720,jmdata=361)
c
      integer n
      integer kpds5u,kpds5v,kpds5t,kpds5z,kpds5q,kpds5rh,kpds5p
      integer kpds5o3,kpds5cl,kpds5pr
      integer kpds5slp
#ifdef DBG
      print *,'lugb,fn=',lugb,fn
      print *,'idim,jdim=',idim,jdim
      print *,'iy,im,id,ih,fh=',iy,im,id,ih,fh
#endif
c
      ijdim=idim*jdim
c
      if(fn(1:4).eq.'    ') then
         print *,'pressure grib file name empty.'
         call abort
      endif
c
      allocate (kpd(25,kdimprs,nfld))
c
      ijmdim=imdata*jmdata
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
c  n=6: o3 (kpds5=154)
c  n=7: cl (kpds5=153)
c  n=8: pr (kpds5=152)
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
         kpds5o3=154
         kpds5cl=153
         kpds5pr=152
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
          elseif (ihyb.eq.1) then
            kpd(6,k,n)=109
          elseif (ihyb.eq.2) then
            kpd(6,k,n)=107
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
        if (ntotal_.gt.1) then
          if (izrtr.eq.0) then
             kpd(5,k,6)=kpds5o3
             kpd(5,k,7)=kpds5cl
             kpd(5,k,8)=kpds5pr
          else
             kpd(5,k,6)=kpds5t     !! temporary
             kpd(5,k,7)=kpds5t
             kpd(5,k,8)=kpds5t
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
         print*,'in read grib n=',n
        do k=1,kdimprs
 82       continue
          print *,'== start reading ',kpd(5,k,n),kpd(6,k,n),kpd(7,k,n)
          call fixrd2(lugb,fn,kpd(1,k,n),
     &               iy,im,id,ih,fh,
     &               imdata,jmdata,
     &               imaxgrb,jmaxgrb,
     &               kpds5,kgds1,
#ifdef G2G
     &               proj,orient,truth,cotru,
     &               delx,dely,rlat1,rlat2,rlon1,rlon2,
     &               dlon,dlat,wlon,rnlat,
#else
     &               rproj,rorient,rtruth,cotru,
     &               rdelx,rdely,rlat1,rlat2,rlon1,rlon2,
#endif
     &               data,iret)
          print*,'finished fixrd2'
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
          if ((inst.eq.1).and.(n.eq.5)) then
             if (isfcp.eq.0) then
                do ij=1,imaxgrb*jmaxgrb
                   data(ij)=data(ij)/9.8
                enddo
             else
                if (k.eq.2) then
                do ij=1,imaxgrb*jmaxgrb
                   data(ij)=data(ij)/9.8
                enddo
                endif
             endif
          endif
c
** when slp is used (isfcp=2), topography is zero.
          if ((isfcp.eq.2).and.(n.eq.5).and.(k.eq.2)) then
            do ij=1,imaxgrb*jmaxgrb
              data(ij)=0.
            enddo
          endif
** when izrtr=1, all tracer is zero.
          if ((izrtr.eq.1).and.(n.ge.6)) then
            do ij=1,imaxgrb*jmaxgrb
              data(ij)=0.
            enddo
          endif
#ifdef G2G
#ifdef DBG
          print *,'imaxgrb,jmaxgrb=',imaxgrb,jmaxgrb
          print *,'kpds5,kgds1=',kpds5,kgds1
          print *,'proj,orient,truth,cotru=',
     1             proj,orient,truth,cotru
          print *,'delx,dely,rlat1,rlat2,rlon1,rlon2=',
     &             delx,dely,rlat1,rlat2,rlon1,rlon2
          print *,'dlon,dlat,wlon,rnlat=',
     1             dlon,dlat,wlon,rnlat
          print *,'data=',(data(ij),ij=1,10)
c
c  spatial interpolation 
c
          print *,'imaxgrb,jmaxgrb,dlon,dlat,wlon,rnlat,idim,jdim=',
     1             imaxgrb,jmaxgrb,dlon,dlat,wlon,rnlat,idim,jdim
#endif
          allocate (dummy1(imaxgrb,jmaxgrb))
          do ij=1,imaxgrb*jmaxgrb
            dummy1(ij,1)=1.
          enddo
          do ij=1,idim*jdim
            dummy2(ij,1)=1.
          enddo
c
          kgau=0
          if(kgds1.eq.4) kgau=1
          call la2ga(data,imaxgrb,jmaxgrb,
     &               abs(dlon),abs(dlat),wlon,rnlat,
     &               out(1,k,n),idim,jdim,.false.,dummy1,dummy2,kgau,
     &               iintp)
c     &               0) !! interpolation
c     &               1) !! take closest grid
cc latitude inversion cc !! only for gaussian 
          if ((rlat1.lt.rlat2).and.(kgau.eq.1)) then
            do i=1,idim*jdim
              tmp(i)=out(i,k,n)
            enddo
            do j=1,jdim
              do i=1,idim
                ij=(j-1)*idim+i
                ij2=(jdim-j)*idim+i
                out(ij2,k,n)=tmp(ij)
              enddo
            enddo
          endif
          deallocate (dummy1)
#else
! add by haiqin li
          do ij=1,ijdim
             out(ij,k,n)=data(ij)
          enddo
#endif
          close(lugb)
        enddo
      enddo
#ifdef G2G
c
c  input non lat/lon, non gaussian and output global not allowed.
c
      if(abs(proj).ne.3..and.abs(proj).ne.4.) then
        print *,'input non lat/lon non gaussian and output global',
     1          ' not allowed.'
        call abort
      endif
#endif
c
      deallocate (data)
      deallocate (kpd)
c
#ifdef G2G
#ifdef DBG 
      do n=1,nfld
        do k=1,kdimprs
          print *,'k=',k,' n=',n
          call nntprt(out(1,k,n),idim,jdim)
        enddo
      enddo
#endif 
#endif
      return
      end
