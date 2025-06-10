#define JCAP jcap
#define LATG latg
#define NPES npes
      subroutine fixio(nsfc,fn,sfcftyp,
     1                 labs,iy,im,id,ih,fh,
     2                 sfcfcs,idim,jdim,ioflag)
c
c  read sfc file
c
c  nsfc   : integer unit number for io (input)
c  fn     : input sfc file name
c  sfcftyp: ['osu1','osu2','noa1' 'vic1'] (input)
c  iy     : integer year of initial state (input/output)
c  im     : integer month of intiial state (input/output)
c  id     : integer day of intiial state (input/output)
c  ih     : integer hour of intiial state (input/output)
c  fh     : real forecast hour (input/output)
c  sfcfcs : real output/input array (input/output)
c  idim   : integer x-dimension of sfcfcs. partial dim for mp (input)
c  jdim   : integer y-dimension of sfcfcs. partial dim for mp (input)
c  ioflag : integer  =0 read, =1 write (input)
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
      character*128 fn
      character*4 sfcftyp
      integer idim,jdim
      integer nsfc
      integer ioflag
c
      integer     idate(4)
      character*8 labs(4)
      real sfcfcs(idim*jdim,*)
c
      integer iy,im,id,ih
      real fh
      integer jy,jm,jd,jh  
      real fhx,rjday
      integer  nrecs,maxlev
c
      integer  m,n,nchar
c
      character*8, allocatable :: gvar(:)
      integer,     allocatable :: lev(:)
c
      if(fn(1:6).ne.'linked') then
#ifdef MP
        if( mype.eq.master ) then
#endif
        call numchar(fn,nchar)
        do n=nchar+1,128
          fn(n:n)=' '
        enddo
        close(nsfc)
        open(unit=nsfc,file=fn(1:nchar),form='unformatted',err=811)
        go to 810
  811   continue
        write(6,*) ' fixio: error in opening file ',fn(1:nchar)
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
  810   continue
#ifdef MP
        endif
#endif
      endif
c
      if(ioflag.eq.0) then
#ifdef MP
        if( mype.eq.master ) then
#endif
        rewind nsfc
        read(nsfc) labs
        read(nsfc) fh,idate
#ifdef MP
        endif
#endif
#ifdef MP
#ifndef RSM
        call mpbcastc(labs,32)
        call mpbcastr(fh,1)
        call mpbcasti(idate,4)
#else
        call rmpbcastc(labs,32)
        call rmpbcastr(fh,1)
        call rmpbcasti(idate,4)
#endif
#endif
        iy=idate(4)
        im=idate(2)
        id=idate(3)
        ih=idate(1)
#ifndef NOPRINT
#ifdef MP
        if(mype.eq.master) then
#endif
   99   format(1h ,'fh, idate=',f10.1,2x,4(1x,i4))
        print *,'fixio field read in from unit=',nsfc
        print 99,fh,idate
#ifdef MP
        endif
#endif
#endif
      else
#ifdef MP
        if( mype.eq.master ) then
#endif
        write(nsfc) labs
#ifdef MP
        endif
#endif
#ifdef RESET_FHOUR
        call vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
        if (jy.eq.iy+1.and.jm.eq.1.and.jd.eq.1.and.jh.eq.0 ) then
          idate(4)=jy
          idate(2)=1
          idate(3)=1
          idate(1)=0
          fhx=0.
        else
          idate(4)=iy
          idate(2)=im
          idate(3)=id
          idate(1)=ih
          fhx=fh
        endif
#else
        idate(4)=iy
        idate(2)=im
        idate(3)=id
        idate(1)=ih
        fhx=fh
#endif
#ifdef MP
        if( mype.eq.master ) then
#endif
        write(nsfc) fhx, idate
#ifdef MP
        endif
#endif
#ifdef DBG
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,'fixio field write to unit=',nsfc
        print 99,fhx,idate
#ifdef MP
        endif
#endif
#endif
      endif
c
c  get sfc file properties
c
      call sfcfld(sfcftyp,0,nrecs,lev,gvar,maxlev)
      allocate (lev(nrecs),gvar(nrecs))
      call sfcfld(sfcftyp,1,nrecs,lev,gvar,maxlev)
c
      m=1
      do n=1,nrecs
         call fixrecio(nsfc,gvar(n),lev(n),sfcfcs(1,m),ioflag)
         m=m+lev(n)
      enddo
c
      deallocate (lev,gvar)
c
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifndef NOPRINT
        print *,' fixio completed.'
#endif
#ifdef MP
      endif
#endif
      return
      end
