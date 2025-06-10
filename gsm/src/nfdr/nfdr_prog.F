#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define PI pi
#include <define.h>
      program nfdr_prog
c
c  nfdr forced by model generate (and observed) forcing
c
c  lugb  integer unit number
c  iy .. integer year of initial time
c  im .. integer month of initial time
c  id .. integer day of initial time
c  ih .. integer hour  of initial time
c  fh .. real forecast hour
c
      implicit none
#include <paramodel.h>
#ifndef RSM
#define LONF2F lonf2_
#define LATG2F latg2_
#else
#define LONF2F igrd12_
#define LATG2F jgrd12_
#endif
#ifndef RSM
#define LONF2F lonf2_
#define LATG2F latg2_
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif /* MP */
#else
#define LONF2F igrd12_
#define LATG2F jgrd12_
#ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif

      integer idim,jdim,ijdim
      parameter(idim=LONF2F/2,jdim=LATG2F*2)
      parameter(ijdim=idim*jdim)
#ifdef RSM
      real, allocatable :: flat(:),flon(:)
      real dlamda0
#endif
#include <comnfdr.h>
#ifdef RSM
#include <rscomloc.h>
#endif
      real*4 st2m2pm(ijdim)
      real*4 srh2m2pm(ijdim)
      real*4 sppt24(ijdim)
      real*4 scloud2pm(ijdim)
      real*4 swind2pm(ijdim)
      real*4 sprcpdur(ijdim)
      real*4 stmax24(ijdim)
      real*4 stmin24(ijdim)
      real*4 srhmax(ijdim)
      real*4 srhmin(ijdim)
c
      real rad
      parameter(rad=pi_/180.)
c
c  namelist variables
c
      integer lugb
      integer iy,im,id,ih
      real fh
c
      character*128 condir
c
      integer ij,k,n
c
      data condir/'    '/
      real rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,
     1          rlftgrd,rbtmgrd
      namelist/namnfdr/condir
      namelist/namloc/
     1  rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd
c
      read(5,namnfdr)
      read(5,namloc)
c
#ifdef RSM
      allocate (flat(ijdim),flon(ijdim))
      call setgrd(flat,flon,rdelx,rdely,dlamda0)
      rlat1=flat(1)/rad
      rlon1=flon(1)/rad
      rlat2=flat(idim*2)/rad
      rlon2=flon(idim*2)/rad
      deallocate (flat,flon)
#endif
      open(unit=33,file='./nfdr.forcing.ieee',form='unformatted',
     1     status='old')
      read(33) st2m2pm
      read(33) srh2m2pm
      read(33) scloud2pm
      read(33) swind2pm
      read(33) sppt24
      read(33) sprcpdur
      read(33) stmax24
      read(33) stmin24
      read(33) srhmax
      read(33) srhmin
      do ij=1,ijdim
        t2m2pm(ij,1)=st2m2pm(ij)        ! degree kelvin
        rh2m2pm(ij,1)=srh2m2pm(ij)      ! percent
        cloud2pm(ij,1)=scloud2pm(ij)    ! fraction
        wind2pm(ij,1)=abs(swind2pm(ij)) ! m/s
        ppt24(ij,1)=sppt24(ij)          ! meter unit
        prcpdur(ij,1)=sprcpdur(ij)      ! in second
        tmax24(ij,1)=stmax24(ij)        ! kelvin
        tmin24(ij,1)=stmin24(ij)        ! kelvin
        rhmax(ij,1)=srhmax(ij)          ! percent
        rhmin(ij,1)=srhmin(ij)          ! percent
      enddo
#ifdef RSM
      call shalfi(t2m2pm,1)
      call shalfi(rh2m2pm,1)
      call shalfi(cloud2pm,1)
      call shalfi(wind2pm,1)
      call shalfi(ppt24,1)
      call shalfi(prcpdur,1)
      call shalfi(tmax24,1)
      call shalfi(tmin24,1)
      call shalfi(rhmax,1)
      call shalfi(rhmin,1)
#else
      call row1ns(t2m2pm)
      call row1ns(rh2m2pm)
      call row1ns(cloud2pm)
      call row1ns(wind2pm)
      call row1ns(ppt24)
      call row1ns(prcpdur)
      call row1ns(tmax24)
      call row1ns(tmin24)
      call row1ns(rhmax)
      call row1ns(rhmin)
#endif
      if(fh.ne.0.) then
        call nfdr_pred(iy,im,id,ih,fh,condir)
      endif
c
      return
      end
