#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define NSOIL nsoil
#define MSUB msub
c
      program sfc0
c
c conversion of various types of sfc files, then
c reads grib sfc files and generates initial sfc file
c
      implicit none
#include <paramodel.h>
c
#ifdef RSM
#define ILOTS igrd12_
#define JLOTS jgrd12_
#else
#define ILOTS lonf2_
#define JLOTS latg2_
#endif
      integer idim,jdim,lsoil,ijdim
      parameter(idim=ILOTS/2,jdim=JLOTS*2,lsoil=lsoil_)
      parameter(ijdim=idim*jdim)
c
#ifdef VICLSM1
      integer nsoil
      parameter (nsoil=nsoil_)
#endif
#ifdef VICLSM2
      integer nsoil,msub,kslmb,nslmb
      parameter (nsoil=nsoil_,msub=msub_)
      parameter (kslmb=lsoil*msub, nslmb=nsoil*msub)
#endif
c
#include <varsfc.h>
#include <comsfcmrg.h>
c
      integer nsfc
      character*128 fnin,fnout
      character*4   sfcftypin
      character*4   sfcftyp
c
      real,        allocatable :: sfcfcsin(:,:)     
      character*8, allocatable :: svar(:)     
      integer,     allocatable :: lev(:)
c
      integer numsfcvin,maxlev,numchi
      integer numsfcsin,j,n
      integer ifcst
c
      character*8 labs(4)
      integer iy,im,id,ih,jy,jm,jd,jh
      real fh
#ifdef ADJSOILW
      real dummy(ILOTS,JLOTS)
      integer idummy,i
#endif
c
#ifdef RSM
c
#include <rscomloc.h>
c
      namelist /namloc/
     1 rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd
     2,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
#endif
c
c
c  namelist variables
c    nsfc     : integer   .. unit number for i/o
c    fnin     : character .. input sfc file name
c    fnout    : character .. output sfc file name
c    sfcftypin: character .. input sfc file type [osu1,osu2,noa1,vic1]
c
      data nsfc/19/
      namelist/namsfc0/nsfc,fnin,fnout,sfcftypin
c
      read(5,namsfc0)
#ifdef RSM
      read(5,namloc)
#endif
c
#ifdef OSULSM1
      sfcftyp='osu1'
#endif
#ifdef OSULSM2
      sfcftyp='osu2'
#endif
#ifdef NOALSM1
      sfcftyp='noa1'
#endif
#ifdef VICLSM1
      sfcftyp='vic1'
#endif
#ifdef VICLSM2
      sfcftyp='vic2'
#endif
c
c  find this file type's record structure and allocate space
c  for reading variables.
c
c    numsfcv is number of surface file variables
c    numsfcs is number of surface file 2-D recs.
c
      call sfcfld(sfcftypin,0,numsfcvin,lev,svar,maxlev)
c
      allocate (lev (numsfcvin))
      allocate (svar(numsfcvin))
c
      call sfcfld(sfcftypin,1,numsfcvin,lev,svar,maxlev)
c
      numsfcsin=0
      do n=1,numsfcvin
        numsfcsin=numsfcsin+lev(n)
      enddo
      allocate (sfcfcsin(ijdim,numsfcsin))
c
      print *,'NUMsfcVARSin=',numsfcvin
      print *,'NUMsfcRECSin=',numsfcsin
c
c  read input sfc file. the file type is sfcftypin and may
c  not be the same as the sfc file type of the model you are running.
c
c
      call fixio(nsfc,fnin,sfcftypin,
     1           labs,iy,im,id,ih,fh,
     2           sfcfcsin,idim,jdim,0)
c
#ifdef GDAS
      call incdte(iy,im,id,ih,
     1            jy,jm,jd,jh,nint(fh))
      iy=jy
      im=jm
      id=jd
      ih=jh
      fh=0.
#else
      if(fh.ne.0.) then
        print *,'fh.ne.0. reset to 0.'
        fh=0.
      endif
#endif
c
c  convert input sfc file type to current model surface type
c
      if(sfcftypin.ne.sfcftyp) then 
        call sfcfcnv(sfcftypin,sfcfcsin,ijdim,numsfcsin)
      else
         do n=1,numsfcsin
            do j=1,ijdim
               sfcfcs(j,n)=sfcfcsin(j,n)
            enddo
         enddo
      endif
c
c
c
c  create initial states by reading grib files and
c  merging it to the sfc file. output is sfcfcs in common
c
      ifcst=0
c
#ifdef ADJSOILW
      do j=1,JLOTS
        do i=1,ILOTS
           dummy(i,j)=0.
        enddo
      enddo
      idummy=24
      call sfc(nsfc,iy,im,id,ih,fh,ifcst,idummy,dummy,dummy)
#else
      call sfc(nsfc,iy,im,id,ih,fh,ifcst)
#endif
c
c
c
c
#ifdef VICLSM1
c
c  initialize some specific parameters and variables used by VIC1
c  one-tile
c
      print *,' in sfc0 before vic1ini'
c
      if(sfcftypin.ne.sfcftyp) then
         call vic1ini(ijdim,im)
      end if
c
#endif
c
#ifdef VICLSM2
c
c  initialize some specific parameters and variables used by VIC2
c  multi-tile
c
      print *,' in sfc0 before vic2ini'
c
      if(sfcftypin.ne.sfcftyp) then
         call vic2ini(ijdim,im)
      end if
c
#endif
c
c  write out
c
c
      call fixio(nsfc,fnout,sfcftyp,
     1           labs,iy,im,id,ih,fh,
     2           sfcfcs,idim,jdim,1)
c
      stop
      end
