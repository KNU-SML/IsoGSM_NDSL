#define JCAP jcap
#define LATG latg
#define NPES npes
      subroutine rdgrbsfc(lugb,idim,jdim,
     &                    numsfcs,numsfcv,ksfc,
     &                    is2g,
     &                    iy,im,id,ih,fh,
     &                    fnmask,fnorog,fnmskg,fn,numgrbs,
     &                    orog,slmask,
     &                    grbfld,lsf)
c
c  read surface grib files
c
c    lugb:  integer unit number 
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
c
      integer lugb,idim,jdim,numsfcs,numsfcv
      integer iy,im,id,ih
      real    fh
      integer numgrbs
      real orog(idim*jdim),slmask(idim*jdim)
      real          grbfld(idim*jdim,numsfcs)
      character*128 fnmask,fnorog,fnmskg
      character*128 fn(numgrbs)
      integer       is2g(numsfcs)
c
      logical       lsf(numsfcs)
      integer       ksfc(numsfcs)
c
      integer nsfc,nsfcv,ngrb
#ifdef DBG
      integer nch
#endif
c
      integer ifp
      data ifp/0/
      save ifp
c
c  read fixed fiedls once at first pass and at fh=0.
c
      if(ifp.eq.0) then
c
c  get land-sea mask on model from mtn program binary output.
c
#ifdef MP
        if(mype.eq.master) then
#endif
        call rdmask(lugb,fnmask,idim*jdim,slmask)
#ifdef MP
        endif
#endif
c
c  get orography on model grid from mtn program binary output.
c
#ifdef MP
        if(mype.eq.master) then
#endif
        call rdorog(lugb,fnorog,idim*jdim,slmask,orog)
#ifdef MP
        endif
#endif
c
        ifp=1
      endif
c
c  read all other grib files
c   note that grbfld is in sfc rec order
c
#ifdef MP
        if(mype.eq.master) then
#endif
      nsfc=1
      do nsfcv=1,numsfcv
        ngrb=is2g(nsfc)
#ifdef DBG
        if(ngrb.ne.9999) then
          call numchar(fn(ngrb),nch)
        print 100,nsfcv,nsfc,fn(ngrb)(1:nch),ngrb,ksfc(nsfc)
  100   format('nsfcv=',i3,' nsfc=',i3,' fn=',a50,' ngrb=',i4,
     1         ' ksfc(nsfc)=',i3)
        endif
#endif
        if(ngrb.ne.9999) then
          if(fn(ngrb)(1:4).ne.'    ') then
            lsf(nsfc)=.true.
            call superead(lugb,fn(ngrb),idim,jdim,slmask,
     1                    ngrb,iy,im,id,ih,fh,
     2                    fnmskg,grbfld(1,nsfc))
          else
            lsf(nsfc)=.false.
          endif
        endif
        nsfc=nsfc+ksfc(nsfc)
      enddo
#ifdef MP
        endif
#endif
c
      return
      end
