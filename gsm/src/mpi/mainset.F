      program mainset
c$$$  main program documentation block
c
c main program: 	mainset
c    prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: this program to get maxima dimension for all partitions
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
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
c     include 'dimset.F'
c     include 'equdiv.F'
c     include 'equdis.F'
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
c
      integer npes,ncol,nrow,n
c
      integer,allocatable::levstr(:),levlen(:)
      integer,allocatable::latstr(:),latlen(:)
      integer,allocatable::lonstr(:),lonlen(:)
      integer,allocatable::lwvstr(:),lwvlen(:)
      integer,allocatable::lntstr(:),lntlen(:)
      integer,allocatable::lnpstr(:),lnplen(:)
      integer,allocatable::lwvdef(:),latdef(:)
c     
      integer jcap,levs,lonf,latg
     *       ,levmax,lonmax,latmax,lntmax,lnpmax,lwvmax
c
      print *,' Enter jcap levs lonf latg npes ncol nrow '
      read(*,*) jcap,levs,lonf,latg,npes,ncol,nrow
c
      allocate(levstr(0:npes-1))
      allocate(latstr(0:npes-1))
      allocate(lonstr(0:npes-1))
      allocate(lwvstr(0:npes-1))
      allocate(lntstr(0:npes-1))
      allocate(lnpstr(0:npes-1))
      allocate(levlen(0:npes-1))
      allocate(latlen(0:npes-1))
      allocate(lonlen(0:npes-1))
      allocate(lwvlen(0:npes-1))
      allocate(lntlen(0:npes-1))
      allocate(lnplen(0:npes-1))
      allocate(lwvdef(jcap+1))
      allocate(latdef(latg/2))
c
      call dimset(jcap,levs,lonf,latg,
     *            npes,ncol,nrow,
     *            levstr,levlen,levmax,
     *            lonstr,lonlen,lonmax,
     *            latstr,latlen,latmax,
     *            lwvstr,lwvlen,lwvmax,
     *            lntstr,lntlen,lntmax,
     *            lnpstr,lnplen,lnpmax,
     *            lwvdef,latdef)
c
      print *,'-------------------- lev ------------------'
      do n=0,npes-1
      write(*,100) n,levstr(n),levlen(n)+levstr(n)-1,levlen(n)
      enddo
      print *,'-------------------- lat ------------------'
      do n=0,npes-1
      write(*,100) n,latstr(n),latlen(n)+latstr(n)-1,latlen(n)
      enddo
      print *,'-------------------- lon ------------------'
      do n=0,npes-1
      write(*,100) n,lonstr(n),lonlen(n)+lonstr(n)-1,lonlen(n)
      enddo
      print *,'-------------------- lwv ------------------'
      do n=0,npes-1
      write(*,100) n,lwvstr(n),lwvlen(n)+lwvstr(n)-1,lwvlen(n)
      enddo
      print *,'-------------------- lnt ------------------'
      do n=0,npes-1
      write(*,100) n,lntstr(n),lntlen(n)+lntstr(n)-1,lntlen(n)
      enddo
      print *,'-------------------- lnp ------------------'
      do n=0,npes-1
      write(*,100) n,lnpstr(n),lnplen(n)+lnpstr(n)-1,lnplen(n)
      enddo
      print *,'-------------------- lwvdef ------------------'
      do n=1,jcap+1
      write(*,100) n,lwvdef(n)
      enddo
      print *,'-------------------- latdef ------------------'
      do n=1,latg/2
      write(*,100) n,latdef(n)
      enddo
100   format(4i10)
c
      print *,' jcap levs lonf latg npes ncol nrow '
      print *,  jcap,levs,lonf,latg,npes,ncol,nrow
      print *,' max of levsp =',levmax
      print *,' max of lonfp =',lonmax
      print *,' max of latgp =',latmax
      print *,' max of jcapp  =',lwvmax
      print *,' max of lntp  =',lntmax
      print *,' max of llnp  =',lnpmax
c
      stop
      end
