#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NSOIL nsoil
#define PI pi
#include <define.h>
      subroutine sfcmrg(grbfld,idim,jdim,iy,im,id,ih,fh,
     &                  orog,slmask,glacier,amxice,
#ifdef ADJSOILW
     &                  intobsp,mdlpcp,mdlrnf,
#endif
     &                  lsf,is2g)
c
c
c  this program takes grib climatology and analysis fields,
c  merge them with the surface varible array in the forecast program.
c
c      implicit none
#include <paramodel.h>
#include <varsfc.h>
#include <vargrb.h>
#include <comsfcmrg.h>
#ifndef RSM
#include <comfgrid.h>
#endif
#ifdef MP
#ifdef RSM
#include <rscommap.h>
#else
#include <comfphys.h>
#endif
#endif
c
c
      integer idim,jdim,iy,im,id,ih
      real    fh
      real    grbfld(idim*jdim,numsfcs)
      real    orog(idim*jdim),slmask(idim*jdim)
      real    glacier(idim*jdim),amxice(idim*jdim)
      logical lsf(numsfcs)
      integer is2g(*)
      character*8 svar
#ifndef MP
      dimension rlat(idim,jdim),rlon(idim,jdim)
#endif
      real rdelx,rdely,dlamda0
      integer i,j
#ifdef USGS
      integer jy,jm,jd,jh
#endif
c
      integer, allocatable :: loismsk(:,:)
#ifdef DBG
      real, allocatable :: rloismsk(:,:)
#endif
c
      real rlapse
      parameter(rlapse=0.65e-2)
c
      real rjday
      integer ij,nsfc,nsfcv,ngrb,k
      character*4 sfcftyp
c
      logical lmerge(numsfcs)
#ifdef ADJSOILW
      real mdlpcp(idim,jdim),mdlrnf(idim,jdim)
#endif
c
c  temporary, set loismsk to land sea mask
c
      allocate (loismsk(idim,jdim))
      do ij=1,idim*jdim
        loismsk(ij,1)=nint(slmask(ij))
      enddo
c
c  quality control of sea ice
c
      call qcsice(grbfld(1,jsli),glacier,amxice,slmask,idim,jdim)
c
c  merge sea ice
c
c    always use slmask as land sea mask, discarding the mask
c    on original sfc file
c
      if (lflgrb(3,iais)) then

      do ij=1,idim*jdim
        sfcfcs(ij,jsli)=slmask(ij)
      enddo
c
c    then add new sea ice
c
      do ij=1,idim*jdim
        if(grbfld(ij,jsli).eq.1.) sfcfcs(ij,jsli)=2.
      enddo
!
      endif
c
c  temporary, set loismsk to land sea sea-ice mask
c
      do ij=1,idim*jdim
        if(sfcfcs(ij,jsli).eq.2.) then
          loismsk(ij,1)=2
        endif
      enddo
c
c  process snow
c   if snow depth available (snow cover analysis not used)
c
      ngrb=is2g(jsno)
      if(lsf(jsno).and.ngrb.ne.9999) then
        ij=1
        do while (grbfld(ij,jsno).ne.-999..and.ij.lt.idim*jdim)
           ij=ij+1
        enddo
        if(ij.eq.idim*jdim) then
          call mrgfld(sfcfcs(1,jsno),grbfld(1,jsno),
     1                idim*jdim,ksfc(jsno),
     2                lsf(jsno),lflgrb(1,ngrb),loismsk)
#ifdef NOALSM1
          do ij=1,idim*jdim
            sfcfcs(ij,jsnd)=sfcfcs(ij,jsno)*5.
          enddo
#endif
        else
c
c  merge snow
c
          call mrgsnow(grbfld(1,jsno),sfcfcs(1,jsno),sfcfcs(1,jtsf),
#ifdef NOALSM1
     1                 sfcfcs(1,jsnd),
#endif
     2                 idim*jdim)
        endif
      endif
c
c  generate loismsk
c        0=sea
c        1=snow free land
c        2=snow free sea-ice
c        3=land with snow
c        4=sea ice with snow mask
c
      do ij=1,idim*jdim
        if(sfcfcs(ij,jsno).gt.0.) then
          if(loismsk(ij,1).eq.1) then
             loismsk(ij,1)=3
          elseif(loismsk(ij,1).eq.2) then
             loismsk(ij,1)=4
          endif
        endif
      enddo
#ifdef DBG
c
      allocate (rloismsk(idim,jdim))
      do ij=1, idim*jdim
        rloismsk(ij,1)=loismsk(ij,1)
      enddo
#ifndef RSM
      call rowsep(rloismsk)
#else
      call shalfo(rloismsk,1)
#endif
      call nntprt(rloismsk,idim,jdim)
#endif
c
c make sst and tg3 correction before merge
c
      if(lsf(jtg3)) then
        call tsfcor(grbfld(1,jtg3),orog,slmask,1.,idim*jdim)
      endif
      if(lsf(jtsf)) then
        call tsfcor(grbfld(1,jtsf),orog,slmask,0.,idim*jdim)
      endif
c
c  now ready to merge 
c
      nsfc=1
      do nsfcv=1,numsfcv
        lmerge(nsfc)=.true.
        if(nsfc.eq.jsli) lmerge(nsfc)=.false.
        if(nsfc.eq.jsno) lmerge(nsfc)=.false.
        if(nsfc.eq.jobp) lmerge(nsfc)=.false.
        nsfc=nsfc+ksfc(nsfc)
      enddo
c
!

      nsfc=1
      do nsfcv=1,numsfcv
        if(lmerge(nsfc)) then      !F: skp merge done above
          ngrb=is2g(nsfc)
!
#ifdef DBG
       if (ngrb.ne.9999) then
            print*,'hoon:',gvar(ngrb),lsf(nsfc),nsfc,ngrb,lflgrb(:,ngrb)
          else
            print*,'hoon:',lsf(nsfc),nsfc,ngrb,lflgrb(:,ngrb)
       endif
#endif
!

          if(lsf(nsfc).and.ngrb.ne.9999) then
            call mrgfld(sfcfcs(1,nsfc),grbfld(1,nsfc),
     1                  idim*jdim,ksfc(nsfc),
     2                  lsf(nsfc),lflgrb(1,ngrb),loismsk)
          endif
        endif
        nsfc=nsfc+ksfc(nsfc)
      enddo
!
#if defined NOA || defined OSULSM2
c
c  quality control of vegataion and soil type
c
      call landtyp(sfcfcs(1,jvet),sfcfcs(1,jsot),sfcfcs(1,jsli),
     &             idim,jdim)
#else
c
c  check vic land surface parameters
c
      call viclandck(idim,jdim)
#endif
#ifdef USGS
c
c  in case of USGS, zor and albedo are computed from vegetation type
c
      call rough(sfcfcs(1,jvet),idim*jdim,sfcfcs(1,jzor))
      call vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
#ifndef MP
#ifdef RSM
      call setgrd(rlat,rlon,rdelx,rdely,dlamda0)
#else
      call lonlat(rlon,rlat,colrad,lonf_,latg_)
      call row1ns(rlon)
      call row1ns(rlat)
      do j=1,latg2_
        do i=1,lonf2_
          rlon(i,j)=rlon(i,j)*180./pi_
          rlat(i,j)=rlat(i,j)*180./pi_
        enddo
      enddo
#endif
#endif
#define RLAT rlat
#ifndef RSM
#ifndef MP
#define RLAT xlat
#endif
#endif
      call setalbedo(sfcfcs(1,jvet),RLAT,idim*jdim,jm,
     1               sfcfcs(1,jalb),sfcfcs(1,jalf))
#endif
c-----------------------------------------------------------
c  quality control by checking max/min over various surfaces
c
      nsfc=1
      do nsfcv=1,numsfcv
        write(svar,100) nsfcv
  100   format('nv=',i3)
        do k=1,ksfc(nsfc)
          call qcmxmn(sfcfcs(1,nsfc+k-1),idim*jdim,loismsk,
     1                vmaxmin(1,1,nsfc),svar)
        enddo
        nsfc=nsfc+ksfc(nsfc)
      enddo

c-----------------------------------------------------------
c
#ifdef NOALSM1
c
c  check snodpth to make sure that it is nonzero if snow is nonzero
c
      call checksndph(sfcfcs(1,jsno),sfcfcs(1,jsnd),idim*jdim)
#endif
#ifdef ADJSOILW
      call adjsoilw(idim,jdim,lsoil_,intobsp,iy,im,id,ih,fh,
     1              grbfld(1,jobp),mdlpcp,
     2              sfcfcs(1,jsli),sfcfcs(1,jtsf),sfcfcs(1,jsno),
     3              sfcfcs(1,jsmc),mdlrnf)
      
#endif
      deallocate (loismsk)
c
#ifdef DBG
      deallocate (rloismsk)
#endif
      return
      end
