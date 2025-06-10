      subroutine dimset(jcap,levs,lonf,latg,
     *                  npes,ncol,nrow,
     *                  levstr,levlen,levmax,
     *                  lonstr,lonlen,lonmax,
     *                  latstr,latlen,latmax,
     *                  lwvstr,lwvlen,lwvmax,
     *                  lntstr,lntlen,lntmax,
     *                  lnpstr,lnplen,lnpmax,
     *                  lwvdef,latdef)
c$$$  subprogram documentation block
c
c subprogram: 	dimset
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: preset all starting point and length for 
c           all pe for global spectral model.
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call dimset(jcap,levs,lonf,latg,
c    *                  npes,ncol,nrow,
c    *                  levstr,levlen,levmax,
c    *                  lonstr,lonlen,lonmax,
c    *                  latstr,latlen,latmax,
c    *                  lwvstr,lwvlen,lwvmax,
c    *                  lntstr,lntlen,lntmax,
c    *                  lnpstr,lnplen,lnpmax,
c    *                  lwvdef,latdef)
c
c    input argument lists:
c	jcap	- integer spectral wavenumber
c	levs	- integer vertical layer number
c	lonf	- integer gaussian grid for longitude
c	latg	- integer gaussian grid for latitude
c	npes	- integer number of pe used: npes=ncol*nrow
c	ncol	- integer number of column
c	nrow	- integer number of nrow
c
c    output argument list:
c	lev*	- integer (npes) related to layers for each pe
c	lon*	- integer (npes) related to longitude for each pe
c	lat*	- integer (npes) related to latitude for each pe
c	lnt*	- integer (npes) related to npes cut of spectral
c	lnp*	- integer (npes) related to nrow cut of spectral
c	lwv*	- integer (npes) related to group of spectral in l
c	*str	- integer (npes) related to each starting point
c	*len	- integer (npes) related to each length
c	*max	- integer related to maximal length of the kind
c	lwvdef	- integer (jcap+1) index of l wave after distribution
c	latdef	- integer (latg/2) index of latitude after distribution
c 
c subprograms called:
c   equdiv	- to compute about equal number of subgroup by division
c   equdis	- to compute about equal number of subgroup by distribution
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer jcap,levs,lonf,latg,npes,ncol,nrow
     *       ,latg2,lonf2,nr,nc,nn,jcaprm,nremain
     *       ,lnp,lx,lh,n
     *       ,levmax,lonmax,latmax,lwvmax,lntmax,lnpmax
     *       ,levpnt,lonpnt,latpnt,lwvpnt,lntpnt,lnppnt
     *       ,levstr,lonstr,latstr,lwvstr,lntstr,lnpstr
     *       ,levlen,lonlen,latlen,lwvlen,lntlen,lnplen
     *       ,lwvdef,latdef
      dimension
     *        levstr(0:npes-1),levlen(0:npes-1)
     *       ,lonstr(0:npes-1),lonlen(0:npes-1)
     *       ,latstr(0:npes-1),latlen(0:npes-1)
     *       ,lwvstr(0:npes-1),lwvlen(0:npes-1)
     *       ,lntstr(0:npes-1),lntlen(0:npes-1)
     *       ,lnpstr(0:npes-1),lnplen(0:npes-1)
     *       ,lwvdef(jcap+1),latdef(latg/2)
      integer, allocatable :: 
     1          levdis(:),londis(:),lntdis(:)
      integer, allocatable :: 
     1          lwvdis(:),latdis(:),lnpdis(:)
c
      allocate (levdis(ncol))
      allocate (londis(ncol))
      allocate (lntdis(ncol))
      allocate (lwvdis(nrow))
      allocate (lnpdis(nrow))
      allocate (latdis(nrow))
c
      call equdiv(levs ,ncol,levdis)
      call equdis( 1,jcap+1,nrow,lwvdis,lwvdef)
c
      lh=0
      do nr=1,nrow
        lnpdis(nr)=0
        lx=lh+1
        lh=lh+lwvdis(nr)
        do n=lx,lh
          lnpdis(nr)=lnpdis(nr)+jcap+2-lwvdef(n)
        enddo
      enddo
c
      latg2=latg/2
      call equdiv(lonf ,ncol,londis)
      call equdis(-1,latg2,nrow,latdis,latdef)
c
      levmax=0
      lonmax=0
      latmax=0
      lwvmax=0
      lntmax=0
      lnpmax=0
c
      latpnt=1
      lwvpnt=1
      lntpnt=1
      lnppnt=1
      n=0
c
      do nr=1,nrow
c
        levpnt=1
        lonpnt=1
        call equdiv(lnpdis(nr),ncol,lntdis)
c
        do nc=1,ncol
c
          levstr(n)=levpnt
          levlen(n)=levdis(nc)
          levpnt=levpnt+levdis(nc)
          levmax=max(levmax,levlen(n))
c
          lonstr(n)=lonpnt
          lonlen(n)=londis(nc)
          lonpnt=lonpnt+londis(nc)
          lonmax=max(lonmax,lonlen(n))
c
          lntstr(n)=lntpnt
          lntlen(n)=lntdis(nc)
          lntpnt=lntpnt+lntdis(nc)
          lntmax=max(lntmax,lntlen(n))
c
          latstr(n)=latpnt
          latlen(n)=latdis(nr)
          latmax=max(latmax,latlen(n))
c
          lwvstr(n)=lwvpnt
          lwvlen(n)=lwvdis(nr)
          lwvmax=max(lwvmax,lwvlen(n))
c
          lnpstr(n)=lnppnt
          lnplen(n)=lnpdis(nr)
          lnpmax=max(lnpmax,lnplen(n))
c
          n=n+1
c
        enddo
c
        if( nr.lt.nrow ) then
          latpnt=latpnt+latdis(nr)
          lwvpnt=lwvpnt+lwvdis(nr)
          lnppnt=lnppnt+lnpdis(nr)
        endif
c
      enddo
c
      do n=0,npes-1
        lwvstr(n)=lwvstr(n)-1
        lnpstr(n)=lnpstr(n)-1
        lntstr(n)=lntstr(n)-1
      enddo
      do n=1,jcap+1
        lwvdef(n)=lwvdef(n)-1
      enddo
c
      levmax=levmax
      lonmax=lonmax
      latmax=latmax*2
      lwvmax=lwvmax	! consider as jcapp
      lntmax=lntmax
      lnpmax=lnpmax
c
      deallocate (levdis)
      deallocate (londis)
      deallocate (lntdis)
      deallocate (lwvdis)
      deallocate (lnpdis)
      deallocate (latdis)
c
      return
      end
