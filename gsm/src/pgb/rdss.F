#define IO io
#define JO jo
#define LEVS levs
#include <define.h>
      subroutine rdss(nss,nflds,ntrace,ncloud,nc,nctop,jcap,
     &           sl,clat,slat,wlat,trig,ifax,eps,epstop,ss,sstop)
#include <paramodel.h>
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)
      parameter(levs=levs_)
c$$$  subprogram documentation block
c
c subprogram:    rdss        read data from a sigma spectral file
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: reads the records of orography, surface pressure,
c           divergence and vorticity, temperature and humidity
c           from a sigma spectral file.  it is assumed that the first
c           two header records of the file have already been read.
c           the gradients of orography and log surface pressure
c           and the wind components are also computed in spectral space.
c           the geopotential of the pressure gradient is computed too.
c           also, some spectral transform utility fields are computed.
c           subprogram trss should be used to transform to grid
c           as well as compute dry temperature and surface pressure
c           and winds and gradients without a cosine latitude factor.
c
c program history log:
c   91-10-31  mark iredell
c
c usage:    call rdss(nss,jcap,nc,nctop,johf,io2,levs,sl,
c   &                 clat,slat,trig,ifax,eps,epstop,ss,sstop)
c
c   input argument list:
c     nss      - integer unit from which to read file
c     jcap     - integer spectral truncation
c     nc       - integer number of spectral coefficients
c     nctop    - integer number of spectral coefficients over top
c     johf    - integer number of latitude pairs in gaussian grid
c     io2    - integer number of valid data points per latitude pair
c     levs     - integer number of levels
c     sl       - real (levs) sigma full level values
c
c   output argument list:
c     clat     - real (johf) cosines of latitude
c     slat     - real (johf) sines of latitude
c     trig     - real (io2) trigonometric quantities for the fft
c     ifax     - integer (20) factors for the fft
c     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))
c     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top
c     ss       - real (nc,6*levs+6) spectral coefs
c     sstop    - real (nctop,6*levs+6) spectral coefs over top
c                (:,1)                             orography
c                (:,2)                             d(orog)/dx
c                (:,3)                             d(orog)/dy
c                (:,4)                             surface pressure
c                (:,5)                             d(lnps)/dx
c                (:,6)                             d(lnps)/dy
c                (:,7:levs+6)                      divergence
c                (:,levs+7:2*levs+6)               vorticity
c                (:,2*levs+7:3*levs+6)             zonal wind
c                (:,3*levs+7:4*levs+6)             meridional wind
c                (:,4*levs+7:5*levs+6)             temperature
c                (:,5*levs+7:5*levs+ntrace*levs+6)  tracers
c                (:,5*levs+7:5*levs+ntrace*levs+6+ncloud*levs)  clouds
c
c   input files:
c     nss      - sigma spectral file
c
c subprograms called:
c   elat         compute latitudes
c   fftfax       compute utility fields for fft
c   gspc         compute utility fields for spectral transform
c   gradq        compute gradient in spectral space
c   dz2uv        compute vector components in spectral space
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension sl(levs),clat(johf),slat(johf),trig(io2),ifax(20)
      real eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)
      real ss(nc,nflds)
      real sstop(nctop,nflds)
c
      real wlat(johf)
      real enn1((jcap+1)*(jcap+2)),elonn1((jcap+1)*(jcap+2)/2)
      real eon((jcap+1)*(jcap+2)/2),eontop(jcap+1)
      real pln((jcap+1)*(jcap+2)/2),plntop(jcap+1)
      real plndx((jcap+1)*(jcap+2)/2),plndy((jcap+1)*(jcap+2)/2)
      real f(io2/2+3,2,3),wfft(io2,2*3)
c
      parameter(g= 9.8000e+0 ,rd= 2.8705e+2 )
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef SCM
      do nf = 1, nflds
        do i = 1, nc
          ss(i,nf) = 0.0
        enddo
      enddo
#else
c  compute utility fields
      call elat(johf,slat,clat,wlat)
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
      call fftfax (io2/2,ifax,trig)
#endif
#ifdef DEFAULT
      call    fax (ifax,io2/2,3)
      call fftrig (trig,io2/2,3)
#endif
      call gspc(jcap,eps,epstop,enn1,elonn1,eon,eontop)
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read sigma spectral data
      nr=(jcap+1)*(jcap+2)
c     print *,'reading sigma record 1'
c     call flush(6)
      read(nss) (ss(i,1),i=1,nr)
c     print *,'reading sigma record 2'
c     call flush(6)
      read(nss) (ss(i,4),i=1,nr)
c     print *,'reading sigma multilevel record 1'
c     call flush(6)
      do k=1,levs
        read(nss) (ss(i,4*levs+6+k),i=1,nr)
      enddo
c     print *,'reading sigma multilevel record 2'
c     call flush(6)
      do k=1,levs
        read(nss) (ss(i,6+k),i=1,nr)
        read(nss) (ss(i,levs+6+k),i=1,nr)
      enddo
c     print *,'reading sigma multilevel record 3'
c     call flush(6)
      do k=1,ntrace*levs+ncloud*levs
        read(nss) (ss(i,5*levs+6+k),i=1,nr)
        print*,' k ss ',k,ss(1,5*levs+6+k)
      enddo
#ifdef SCM
      do k=1,levs
        read(nss) (ss(i,2*levs+6+k),i=1,nr)  ! work space for SCM
cyh        print*,'Vvel',k,ss(1,2*levs+6+k)
      enddo
#endif
#ifndef SCM
      do k=1,nflds
        do l=0,jcap
          sstop(2*l+1,k)=0.
          sstop(2*l+2,k)=0.
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute gradients and winds
      call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,1),
     &           ss(1,2),ss(1,3),sstop(1,3))
      call gradq(jcap,enn1,elonn1,eon,eontop,ss(1,4),
     &           ss(1,5),ss(1,6),sstop(1,6))
      do k=1,levs
        call dz2uv(jcap,enn1,elonn1,eon,eontop,ss(1,6+k),ss(1,levs+6+k),
     &             ss(1,2*levs+6+k),ss(1,3*levs+6+k),
     &             sstop(1,2*levs+6+k),sstop(1,3*levs+6+k))
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      return
      end
