#define IOSGB iosgb
#define JCAP jcap
#define JOSGB josgb
#define LEVS levs
#define RV rv
#define RD rd
#define G g
#define CP cp
#define HVAP hvap
#include <define.h>
      subroutine trss(trig,ifax,eps,epstop,ss,sstop,coslat,sinlat,
     $                mflds,
     1                f,j,si,sl)
c
#include <paramodel.h>
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)
      parameter(jcap=jcap_,levs=levs_)
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)
c
c  number of spectral conversion fields
c
c      parameter(mflds=18*levs_+12)
c
      dimension trig(io2),ifax(20)
      real eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)
      real ss(nc,mflds),sstop(nctop,mflds)
      real f(io22,mflds)
c
      real si(levs+1),sl(levs)
c
      real pln((jcap+1)*(jcap+2)/2),plntop(jcap+1)
      real wfft(io22,2*mflds)
      parameter(fv_=rv_/rd_-1.)
c
c$$$  subprogram documentation block
c
c subprogram:    trss        transform spectral to grid
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: transforms spectral to gridded data on a latitude pair
c           and computes dry temperature and surface pressure
c           and winds and gradients without a cosine latitude factor.
c           subprogram rdss should be called already
c           to read spectral data and initialize utility fields.
c           this subprogram can be called from a multiprocessed segment.
c
c program history log:
c   91-10-31  mark iredell
c   03-09-16  kana
c
c usage:    call trss(trig,ifax,eps,epstop,ss,sstop,coslat,sinlat,f,j)
c
c   input argument list:
c     trig     - real (io2) trigonometric quantities for the fft
c     ifax     - integer (20) factors for the fft
c     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))
c     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top
c     ss       - real (nc,6*levs+6) spectral coefs
c     sstop    - real (nctop,6*levs+6) spectral coefs over top
c     coslat   - real cosine of latitude of the northern latitude
c     sinlat   - real sine of latitude of the northern latitude
c
c   output argument list:
c     f     - real (io22,6*levs+6) gridded data
c            (:, 1:levs              ) vorticity
c            (:, 1*levs+ 1: 2*levs   ) divergence
c            (:, 2*levs+ 1: 3*levs   ) temperature
c            (:, 3*levs+ 1: 4*levs   ) specific humidity
c            (:, 4*levs+ 1: 5*levs   ) u-component of wind
c            (:, 5*levs+ 1: 6*levs   ) v-component of wind
c            (:, 6*levs+ 1)            surface pressure
c            (:, 6*levs+ 2)            surface geopotential
c            (:, 6*levs+ 3)            d(lnps)/dx
c            (:, 6*levs+ 4)            d(lnps)/dy
c            (:, 6*levs+ 5)            d(orog)/dx
c            (:, 6*levs+ 6)            d(orog)/dy
c
c            (:, 6*levs+ 7: 7*levs+ 6) geopotential energy
c            (:, 7*levs+ 7: 8*levs+ 6) x-gradient of tv
c            (:, 8*levs+ 7: 9*levs+ 6) y-gradient of tv
c            (:, 9*levs+ 7:10*levs+ 6) x-gradient of q
c            (:,10*levs+ 7:11*levs+ 6) y-gradient of q
c            (:,11*levs+ 7:12*levs+ 6) x-derivative of gz
c            (:,12*levs+ 7:13*levs+ 6) y-derivative of gz
c
c            (:,13*levs+ 7:14*levs+ 6) enthalpy (Cp*T)
c            (:,14*levs+ 7)            vertical integral of CpT
c            (:,14*levs+ 8:15*levs+ 7) flux divergence of cpt
c            (:,15*levs+ 8)            vertically integrated flux div of CpT
c
c            (:,15*levs+ 9)            vertically integrated q
c            (:,15*levs+10:16*levs+ 9) flux divergence of q
c            (:,16*levs+10)            vertical integral of the flux div of q
c
c            (:,16*levs+11)            vertically integrated geopotential energy
c            (:,16*levs+12:17*levs+11) flux divergence of gz
c            (:,17*levs+12)            vertically integrated flux div of gz
c
c            (:,17*levs+13:18*levs+12) omega
c            (:,18*levs+13:18*levs+ntrace*levs+13) tracers
c            (:,18*levs+13:18*levs+ntrace*levs+13+ncloud*levs) clouds
c
c subprograms called:
c   pleg         compute associated legendre polynomials
c   psynth       synthesize fourier from spectral coefficients 
c   rfftmlt      fast fourier transform
c
c attributes:
c   language: cray fortran
c
c$$$
c  transform spectral coefficients to fourier coefficients
c
      call pleg(jcap,sinlat,coslat,eps,epstop,pln,plntop)
c
      call psynth(jcap,io22/2,nc,nctop,mflds,pln,plntop,ss,sstop,f)
c
c  transform fourier coefficients to gridded data
c
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
        call rfftmlt(f,wfft,trig,ifax,1,io22/2,io2/2,2*mflds,1)
#endif
#ifdef DEFAULT
      call fft99m (f,wfft,trig,ifax,1,io22/2,io2/2,2*mflds,1)
#endif
c
c  move southern hemisphere latitude after northern hemisphere latitude
c
      do k=1,mflds
#ifdef CRAY_THREAD
cdir$ ivdep
#endif
        do i=1,io2/2
          f(io2/2+i,k)=f(io22/2+i,k)
        enddo
        write(6,*) "trss f=",f(1,k),"ss=",ss(1,k),"k=",k
      enddo
c
c  compute dry temperature from virtual temperature
c  and surface pressure from log surface pressure
c  and divide gradients and winds by cosine of latitude.
c
      do k=1,levs
        do i=1,io2
          f(i, 2*levs  +k)=f(i, 2*levs  +k)/(1.+fv_*f(i,3*levs+k))
          f(i,13*levs+6+k)=f(i, 2*levs  +k)*cp_
          f(i, 6*levs+6+k)=f(i, 6*levs+6+k)*g_
          f(i, 4*levs  +k)=f(i, 4*levs+  k)/coslat
          f(i, 5*levs  +k)=f(i, 5*levs+  k)/coslat
          f(i, 7*levs+6+k)=f(i, 7*levs+6+k)/coslat
          f(i, 8*levs+6+k)=f(i, 8*levs+6+k)/coslat
          f(i, 9*levs+6+k)=f(i, 9*levs+6+k)/coslat
          f(i,10*levs+6+k)=f(i,10*levs+6+k)/coslat
          f(i,11*levs+6+k)=f(i,11*levs+6+k)/coslat*g_
          f(i,12*levs+6+k)=f(i,12*levs+6+k)/coslat*g_
        enddo
      enddo
c
      do i=1,io2
        f(i,6*levs+1)=exp(f(i,6*levs+1))
        f(i,6*levs+3)=f(i,6*levs+3)/coslat
        f(i,6*levs+4)=f(i,6*levs+4)/coslat
        f(i,6*levs+5)=f(i,6*levs+5)/coslat
        f(i,6*levs+6)=f(i,6*levs+6)/coslat
      enddo
c
c  compute flux divergence of temperature and humidity
c
      do k=1,levs
        do i=1,io2
c dt/dx
          f(i,7*levs+6+k)=(f(i,7*levs+6+k)-f(i,2*levs+k)*
     1         fv_*f(i,9*levs+6+k))/(1.+fv_*f(i,2*levs+k))*cp_
c dt/dy
          f(i,8*levs+6+k)=(f(i,8*levs+6+k)-f(i,2*levs+k)*
     1         fv_*f(i,10*levs+6+k))/(1.+fv_*f(i,2*levs+k))*cp_
c div(v*t)
          f(i,14*levs+7+k)=f(i,4*levs+k)*f(i,7*levs+6+k)+
     1                     f(i,5*levs+k)*f(i,8*levs+6+k)+
     2                     f(i,13*levs+6+k)*f(i,levs+k)
c div(v*q)
          f(i,15*levs+9+k)=f(i,4*levs+k)*f(i, 9*levs+6+k)+
     1                     f(i,5*levs+k)*f(i,10*levs+6+k)+
     2                     f(i,3*levs+k)*f(i,levs+k)
c div(v*z)
          f(i,16*levs+11+k)=f(i,4*levs+k)*f(i,11*levs+6+k)+
     1                      f(i,5*levs+k)*f(i,12*levs+6+k)+
     2                      f(i,6*levs+6+k)*f(i,levs+k)
        enddo
      enddo
      call omega(io2,io22,levs,si,sl,f(1,6*levs+1),
     1           f(1,6*levs+3),f(1,6*levs+4),
     2           f(1,  levs+1),f(1,4*levs+1),f(1,5*levs+1),
     3           f(1,17*levs+13))
c
      do i=1,io2
         f(i,14*levs+ 7)=0.
         f(i,15*levs+ 9)=0.
         f(i,16*levs+11)=0.
         f(i,15*levs+ 8)=0.
         f(i,16*levs+10)=0.
         f(i,17*levs+12)=0.
      enddo
      do i=1,io2
        do k=1,levs
          f(i,14*levs+ 7)=f(i,14*levs+ 7)+f(i,13*levs+ 6+k)*
     1                    (si(k)-si(k+1))
          f(i,15*levs+ 9)=f(i,15*levs+ 9)+f(i, 3*levs   +k)*
     1                    (si(k)-si(k+1))
          f(i,16*levs+11)=f(i,16*levs+11)+f(i, 6*levs+ 6+k)*
     1                    (si(k)-si(k+1))
          f(i,15*levs+ 8)=f(i,15*levs+ 8)+f(i,14*levs+ 7+k)*
     1                    (si(k)-si(k+1))
          f(i,16*levs+10)=f(i,16*levs+10)+f(i,15*levs+ 9+k)*
     1                    (si(k)-si(k+1))
          f(i,17*levs+12)=f(i,17*levs+12)+f(i,16*levs+11+k)*
     1                    (si(k)-si(k+1))
        enddo
      enddo
      do i=1,io2
         f(i,14*levs+ 7)=f(i,14*levs+ 7)*f(i,6*levs+1)/g_*1.e3
         f(i,15*levs+ 9)=f(i,15*levs+ 9)*f(i,6*levs+1)/g_*1.e3
         f(i,16*levs+11)=f(i,16*levs+11)*f(i,6*levs+1)/g_*1.e3
         f(i,15*levs+ 8)=f(i,15*levs+ 8)*f(i,6*levs+1)/g_*1.e3
         f(i,16*levs+10)=f(i,16*levs+10)*f(i,6*levs+1)/g_*1.e3
         f(i,17*levs+12)=f(i,17*levs+12)*f(i,6*levs+1)/g_*1.e3
      enddo
c
      return
      end
