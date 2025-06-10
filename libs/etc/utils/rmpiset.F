      program rmpiset
c$$$  main program documentation block
c
c main program: 	rmpiset
c    prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: this program to get maxima dimension for all partitions
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c    00-03-09  henry juang 	update with symmetry distribution
c                               prepare for reduce grid.
c
c input file lists:
c	unit 5	- standard input
c
c output file list:
c	unit 6	- standar output
c 
c subprograms called:
c   dimset	- to compute all maximal dimension, starting
c		  point and length of each pe
c
c     include 'rdimset.F'
c     include 'equdiv.F'
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
c
      integer npes,ncol,nrow,n,nc,igrd1,jgrd1,iwav1,jwav1
c
      integer levr,igrd,jgrd
     *       ,levmax,lonmax,latmax,lntmax,lnpmax,lwvmax
c
      real pesx
      data ncol/0/, nrow/0/
c
      namelist /namset/ levr,igrd,jgrd,npes,ncol,nrow
c
      read(5,namset)
c
      igrd1=igrd+1
      jgrd1=jgrd+1
      iwav1=int((igrd-12)/3)*2+1
      jwav1=int((jgrd-12)/3)*2+1
c
      if( ncol.eq.0 ) then
        pesx=npes
        nc=sqrt(pesx)
        do n=nc,2,-1
          if( mod(npes,n).eq.0 ) then
            ncol=n
            go to 1234
          endif
        enddo
        ncol=1
      else
        if( mod(npes,ncol).ne.0 ) then
          print *,' user provided npes=',npes
          print *,' with an invalid ncol=',ncol
          call abort
          stop
        endif
      endif
1234  nrow=npes/ncol   
c
      call rdimset(iwav1,jwav1,levr,igrd1,jgrd1,
     *            npes,ncol,nrow,
     *            levmax,
     *            lonmax,
     *            latmax,
     *            lwvmax,
     *            lntmax,
     *            lnpmax)
c
      write(*,107) levmax
      write(*,108) lonmax
      write(*,109) latmax
      write(*,110) lwvmax
      write(*,111) lntmax
      write(*,112) lnpmax
c
 101  format('#define MP')
 102  format('#define mpi_real mpi_real8')
 103  format('#define _mpi_real_ 8')
 104  format('#define _npes_ ',i4)
 105  format('#define _ncol_ ',i4)
 106  format('#define _nrow_ ',i4)
 107  format('#define _levrp_ ',i3)
 108  format('#define _igrd1p_ ',i4)
 109  format('#define _jgrd1p_ ',i4)
 110  format('#define _iwav1p_ ',i3)
 111  format('#define _lnwavp_ ',i6)
 112  format('#define _llwavp_ ',i6)
 
      stop
      end
      subroutine rdimset(iwav1,jwav1,levr,igrd1,jgrd1,
     *                  npes,ncol,nrow,
     *                  levmax,
     *                  lonmax,
     *                  latmax,
     *                  lwvmax,
     *                  lntmax,
     *                  lnpmax)
c$$$  subprogram documentation block
c subprogram: 	rdimset
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: preset all starting point and length for 
c           all pe for global spectral model.
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call rdimset(levr,igrd1,jgrd1,
c    *                  npes,ncol,nrow,
c    *                  levstr,levlen,levmax,
c    *                  lonstr,lonlen,lonmax,
c    *                  latstr,latlen,latmax,
c    *                  lwvstr,lwvlen,lwvmax,
c    *                  lntstr,lntlen,lntmax,
c    *                  lnpstr,lnplen,lnpmax)
c
c    input argument lists:
c	iwav1	- integer spectral wavenumber
c	levr	- integer vertical layer number
c	igrd1	- integer gaussian grid for longitude
c	jgrd1	- integer gaussian grid for latitude
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
c 
c subprograms called:
c   equdiv	- to compute about equal number of subgroup by division
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer levr,igrd1,jgrd1,npes,ncol,nrow,iwav1,jwav1
     *       ,latg2,lonf2,nr,nc,nn,jcaprm,nremain
     *       ,lnp,n2a,lx,ll,lh,n,npesx
     *       ,levmax,lonmax,latmax,lwvmax,lntmax,lnpmax
     *       ,levpnt,lonpnt,latpnt,lwvpnt,lntpnt,lnppnt
      parameter(npesx=10000)
      integer
     *        levstr(0:npesx-1),levlen(0:npesx-1)
     *       ,lonstr(0:npesx-1),lonlen(0:npesx-1)
     *       ,latstr(0:npesx-1),latlen(0:npesx-1)
     *       ,lwvstr(0:npesx-1),lwvlen(0:npesx-1)
     *       ,lntstr(0:npesx-1),lntlen(0:npesx-1)
     *       ,lnpstr(0:npesx-1),lnplen(0:npesx-1)
      integer
     1          levdis(npesx),londis(npesx),lntdis(npesx)
     1       ,  lwvdis(npesx),latdis(npesx),lnpdis(npesx)
     1       ,  lwvdef(npesx)
c
      call equdiv(levr ,ncol,levdis)
      call equdiv(iwav1 ,nrow,lwvdis)
c
      do nr=1,nrow
        lnpdis(nr)=lwvdis(nr)*jwav1
      enddo                  

c
      latg2=jgrd1/2
      call equdiv(igrd1 ,ncol,londis)
      call equdiv(latg2 ,nrow,latdis)
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
c
      levmax=levmax
      lonmax=lonmax
      latmax=latmax*2
      lwvmax=lwvmax	! consider as iwv1p
      lntmax=lntmax
      lnpmax=lnpmax
c
      return
      end

      subroutine equdiv(len,ncut,lenarr)
c$$$  subprogram documentation block
c
c subprogram: 	equdiv
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: cut len into ncut pieces with load balancing
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	equdiv(len,ncut,lenarr)
c
c    input argument lists:
c	len	- integer total length 
c	ncut	- integer number of subgroup
c
c    output argument list:
c	lenarr	- integer (ncut) length of each subgroup
c 
c subprograms called: none
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer len,ncut,lenarr
     *       ,n0,n1,n
      dimension lenarr(ncut)
      n0=len/ncut
      n1=mod(len,ncut)
      do n=1,n1
        lenarr(n)=n0+1
      enddo
      do n=n1+1,ncut
        lenarr(n)=n0
      enddo
      return
      end

