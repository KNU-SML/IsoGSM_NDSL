#define JCAP jcap
#define LATG latg
#define NPES npes
      subroutine superead(lugb,fn,idim,jdim,slmask,
     1                    indx,iy,im,id,ih,fh,
     2                    fnmskg,out)
c
c  read grib file, interpolates in time and to model grid.
c  complex interpolation by using masks are applied
c  control of interpolation is specified in inlude/vargrb.h
c
      implicit none
#include <paramodel.h>
#ifdef MP
#ifndef RSM
#include <commpi.h>
#else
#include <comrmpi.h>
#endif
#endif
#include <vargrb.h>
c
c  input:
c
c    lugb    int                   .. unit number for reading grib
c    fn      char                  .. grib file name
c    idim    int                   .. output x-dimension (model x-dimension)
c    jdim    int                   .. output y-dimension (model y-dimension)
c    slmask  real                  .. land sea mask on model grid
c    iy      int                   .. 4-digit year
c    im      int                   .. month
c    id      int                   .. day
c    ih      int                   .. hour
c    fh      real                  .. forecast hour
c    indx    int                   .. pointer of grb fields (see vargrb.h)
c    fnmskg character             .. hi-res lat/lon land sea global mask 
c
c  output:
c   
c    out    real array (idim*jdim,kdim) .. output field
c    iclmtyp integer               .. =1 yearly mean, 4=seasonal mean,
c                                     =12 monthly mean, -1=not a climatology
c
c  the map related parameters appear in this program are for grib
c  files and not those of the forecast program
c
      integer       lugb
      character*128 fn
      integer       idim,jdim
      real          slmask(idim*jdim)
      integer       iy,im,id,ih
      real          fh
      integer       indx
      character*128 fnmskg
      real          out(idim*jdim,*)
      logical lmask
c
      real undef
      parameter(undef=1.e30)
c
      real, allocatable :: data(:),rslmsk(:)
c
      integer k,kgau,ij,kpds5,kgds1,ijmdim,imaxgrb,jmaxgrb
      real proj,orient,truth,cotru
      real delx,dely,rlat1,rlat2,rlon1,rlon2
      real dlon,dlat,wlon,rnlat
      integer iret
      real dummy
      real, allocatable :: tmp(:)
      integer i,j,ij2
c
#ifdef DBG
      print *,'lugb,fn=',lugb,fn
      print *,'idim,jdim=',idim,jdim,kgrb(indx)
      print *,'iy,im,id,ih,fh=',iy,im,id,ih,fh
      print *,'indx=',indx
#endif
c
      if(fn(1:4).eq.'    ') then
         print *,gvar(indx),' file name empty.'
         return
      endif
c
      ijmdim=imdata(indx)*jmdata(indx)
      allocate (data(ijmdim))
c
      lmask=mask(indx).gt.0
c
c  read grib file
c
c  convert land mask to ocean mask if mask=2
c
      if(mask(indx).eq.2) then
         do ij=1,idim*jdim
           slmask(ij)=1.-slmask(ij)
         enddo
      endif
c
      do k=1,kgrb(indx)
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,'========== start reading ',gvar(indx),' ========='
        call fixrd(lugb,fn,kpd(1,k,indx),
     &             iy,im,id,ih,fh,
     &             imdata(indx),jmdata(indx),
     &             imaxgrb,jmaxgrb,
     &             kpds5,kgds1,
     &             proj,orient,truth,cotru,
     &             delx,dely,rlat1,rlat2,rlon1,rlon2,
     &             dlon,dlat,wlon,rnlat,
     &             data,iret)
        if(iret.gt.0) then
          print *,gvar(indx),' read failed at k=',k
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
#ifdef DBG
        print *,'imaxgrb,jmaxgrb=',imaxgrb,jmaxgrb
        print *,'kpds5,kgds1=',kpds5,kgds1
        print *,'proj,orient,truth,cotru=',proj,orient,truth,cotru
        print *,'delx,dely,rlat1,rlat2,rlon1,rlon2=',
     &           delx,dely,rlat1,rlat2,rlon1,rlon2
        print *,'dlon,dlat,wlon,rnlat=',dlon,dlat,wlon,rnlat
        print *,'data=',(data(ij),ij=1,10)
        print *,'lmask,limask(indx),mask(indx)=',
     &           lmask,limask(indx),mask(indx)
#endif
        close(lugb)
#ifdef MP
        endif
#endif
cc latitude inversion cc !! only for gaussian
          if ((rlat1.lt.rlat2).and.(kgds1.eq.4)) then
            print*,'inversed!',indx,ijmdim,imdata(indx),jmdata(indx)
            allocate (tmp(ijmdim))
            do i=1,ijmdim
              tmp(i)=data(i)
            enddo
            do j=1,jmaxgrb
              do i=1,imaxgrb
                ij=(j-1)*imaxgrb+i
                ij2=(jmaxgrb-j)*imaxgrb+i
                data(ij2)=tmp(ij)
              enddo
            enddo
            deallocate (tmp)
          endif
c
c  define mask on input grid 
c
c  create slmask over input grid using ready made hires lat/lon mask
c  when no mask info is available from grib file itself
c
        if(k.eq.1) then
          allocate (rslmsk(imaxgrb*jmaxgrb))
        endif
        if(lmask) then
          if(.not.limask(indx)) then
            if(k.eq.1) then
              call ll2msk(lugb,fnmskg,
     &                    0,rslmsk,imaxgrb,jmaxgrb,
     &                    proj,orient,truth,cotru,
     &                    delx,dely,rlat1,rlat2,rlon1,rlon2)
c  convert land mask to ocean mask if mask=2
              if(mask(indx).eq.2) then
                do ij=1,imaxgrb*jmaxgrb
                  rslmsk(ij)=1.-rslmsk(ij)
                enddo
              endif
            endif
          else
            call fld2msk(data,imaxgrb*jmaxgrb,
     1                   cvalin(indx),condin(indx),rslmsk)
c  convert land mask to ocean mask if mask=2
            if(mask(indx).eq.2) then
              do ij=1,imaxgrb*jmaxgrb
                rslmsk(ij)=1.-rslmsk(ij)
              enddo
            endif
          endif
        else
          do ij=1,imaxgrb*jmaxgrb
            rslmsk(ij)=1.
          enddo
          if(mask(indx).eq.2) then
            do ij=1,imaxgrb*jmaxgrb
              rslmsk(ij)=1.-rslmsk(ij)
            enddo
          endif
        endif
c
c  spacial interpolation 
c
#ifdef RSM
        call ll2rsm(data,imaxgrb,jmaxgrb,
     &              inttyp(indx),out(1,k),idim,jdim,
     1              lmask,rslmsk,slmask,
     2              proj,orient,truth,cotru,
     3              delx,dely,rlat1,rlat2,rlon1,rlon2)
#else
c
        kgau=0
        if(kgds1.eq.4) kgau=1
        call la2ga(data,imaxgrb,jmaxgrb,
     &             abs(dlon),abs(dlat),wlon,rnlat,
     &             out(1,k),idim,jdim,lmask,rslmsk,slmask,kgau,
     &             inttyp(indx))
#endif
      enddo
c
      deallocate (data)
      deallocate (rslmsk)
c
c  get original mask value if mask=2
c
      if(mask(indx).eq.2) then
        do ij=1,idim*jdim
          slmask(ij)=1.-slmask(ij)
        enddo
      endif
c
c  scale the value if necessary
c
       if(scale(indx).ne.1.) then
         do k=1,kgrb(indx)
           do ij=1,idim*jdim
             out(ij,k)=out(ij,k)*scale(indx)
           enddo
         enddo
       endif
c
c  process mask output (0 or 1 field)
c 
       if(lomask(indx)) then
         call fld2msk(out,idim*jdim,cvalout(indx),condout(indx),out)
       endif
c
#ifdef DBG 
      do k=1,kgrb(indx)
        print *,gvar(indx),' for k=',k
        call nntprt(out(1,k),idim,jdim)
      enddo
#endif 
      return
      end
